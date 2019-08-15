
# install.packages("remotes")
#remotes::install_github("centerforopenscience/osfr")

library(osfr)

library(utils)
## Packages
require(dplyr)
require(tidyverse)
require(fmri)
require(boot)
require(ggpubr)


#see API token set up at: http://centerforopenscience.github.io/osfr/articles/auth.html
get_data=0

if(get_data==1)
{
  user_renviron = path.expand(file.path("/Volumes/PSYHOME/PSYRES/pthompson/DVMB/fTCD_glm_project", ".Renviron"))
  if(!file.exists(user_renviron)) # check to see if the file already exists
    file.create(user_renviron)
  file.edit(user_renviron) #paste in the PAT token created earlier into file like this: OSF_PAT='insert token token'
  
  osf_auth(token=Sys.getenv("OSF_PAT"))
  
  lb_project <- osf_retrieve_node("hfn2j")
  
  osf_retrieve_file("https://osf.io/5kq42/") %>%
    osf_download() %>% unzip('Chpt4_fTCD_WordGen_rawdata.zip')
}
#########################################################################################
#########################################################################################

fTCD_glm3<-function(path,order)
{
  
  
  #---------------------------------------------------------------------------------------------------------------#
  
  filename1<-list.files(path,pattern = '.exp')
  ## Set parameters
  samplingrate <- 25 # Sampling rate after downsampling. Raw data is 100Hz, we take 1 in every 4 samples
  heartratemax <- 125
  
  
  glm.data<-data.frame(matrix(NA,nrow=length(filename1),ncol=(((order+3)*2)+3)))
  names(glm.data)<-c('ID',paste0('Lparam',(1:(order+3))),'L_HRF',paste0('Rparam',(1:(order+3))),'R_HRF')
  
  ts_summary_data<-as.data.frame(matrix(NA,nrow=length(filename1),ncol=5))
  names(ts_summary_data)<-c("ID","peakL","bmeanL","peakR","bmeanR")
  
  
  for(j in 1:length(filename1))
  {
    print(filename1[j])
    #########################################################################################
    # fTCD preprocessing for GLM analysis
    # 
    # Created by z.woodhead 30th July 2019 
    #########################################################################################
    # This script takes a raw .exp datafile and preprocesses it ready for GLM analysis:
    #   - It creates a box car function showing when the task was ON or OFF
    #   - It normalises the fTCD signal to a mean of 100
    #   - It performs heart beat integration
    #   - It saves the processed data into a .csv file
    
    ## Read in raw data
    
    myfile <- filename1[j]
    mydata<-read.table(myfile, skip = 6,  header =FALSE, sep ='\t')
    
    wantcols = c(2,3,4,7) #sec, L, R,marker #select columns of interest to put in shortdat
    shortdat = data.frame(mydata[,wantcols])
    rawdata = filter(shortdat, row_number() %% 4 == 0) # downsample to 25 Hz by taking every 4th point
    allpts = nrow(rawdata) # total N points in long file
    rawdata[,1] = (seq(from=1,to=allpts*4,by=4)-1)/100 #create 1st column which is time in seconds from start
    colnames(rawdata) = c("sec","L","R","marker")
    
    #----------------------------------------------------------
    ## Find markers; place where 'marker' column goes from low to high value
    
    mylen = nrow(rawdata); # Number of timepoints in filtered data (rawdata)
    markerplus = c(rawdata$marker[1] ,rawdata$marker); # create vectors with offset of one
    markerchan = c(rawdata$marker,0); 
    markersub = markerchan - markerplus; # start of marker indicated by large difference between consecutive data points
    meanmarker <- mean(rawdata$marker) # We will identify big changes in marker value that are > 5 sds
    markersize <- meanmarker+4*sd(rawdata$marker)
    origmarkerlist = which(markersub>markersize)
    norigmarkers = length(origmarkerlist)
    
    # Stimulus timings: Word Gen starts 5 seconds after marker and continues for 20 seconds (including REPORT phase)
    # Edit: stim1 models covert word generation, which starts 5 seconds after marker and continutes for 15 seconds
    # stim2 models overt word reporting, which starts 20 seconds after marker and continues for 5 seconds
    stim1_delay_sec <- 5
    stim1_delay_samples <- stim1_delay_sec * samplingrate
    stim1_length_sec <- 15
    stim1_length_samples <- stim1_length_sec * samplingrate
    
    stim2_delay_sec <- 20
    stim2_delay_samples <- stim2_delay_sec * samplingrate
    stim2_length_sec <- 5
    stim2_length_samples <- stim2_length_sec * samplingrate
    
    rest_length_sec <- 30
    rest_length_samples <- rest_length_sec * samplingrate
    
    rawdata$stim1_on <- 0
    rawdata$stim2_on <- 0
    for (m in 1:norigmarkers){
      rawdata$stim1_on[(origmarkerlist[m]+stim1_delay_samples):(origmarkerlist[m]+stim1_delay_samples+stim1_length_samples)] <- 1
      rawdata$stim2_on[(origmarkerlist[m]+stim2_delay_samples):(origmarkerlist[m]+stim2_delay_samples+stim2_length_samples)] <- 1
    }
    
    #----------------------------------------------------------
    # Data normalisation
    
    meanL=mean(rawdata$L)
    meanR=mean(rawdata$R)
    rawdata$normal_L=rawdata$L/meanL * 100 #last dim of myepoched is 2 for the normalised data
    rawdata$normal_R=rawdata$R/meanR * 100
    
    
    #----------------------------------------------------------
    # Heartbeat integration
    peaklist=numeric(0)
    pdiff=numeric(0)
    badp=numeric(0)
    
    # Look through every sample from 6, to number of samples minus 6
    for(i in seq(6,mylen-6))
    {if(
      (rawdata$L[i] > rawdata$L[i-5])
      & (rawdata$L[i] > rawdata$L[i-4])
      & (rawdata$L[i] > rawdata$L[i-3])
      & (rawdata$L[i] > rawdata$L[i-2])
      & (rawdata$L[i] > rawdata$L[i-1])
      & (rawdata$L[i] > rawdata$L[i+1])
      & (rawdata$L[i] > rawdata$L[i+2])
      & (rawdata$L[i] > rawdata$L[i+3])
      & (rawdata$L[i]> rawdata$L[i+4])
      & (rawdata$L[i]> rawdata$L[i+5]))
    {peaklist=c(peaklist,i)
    }
    }
    
    # Check that the heartbeats are spaced by far enough!
    peakdiffmin = 60/heartratemax * samplingrate
    pdiff <- peaklist[2:length(peaklist)]-peaklist[1:(length(peaklist)-1)] # pdiff is a list of the number of samples between peaks
    badp<-which(pdiff<peakdiffmin) # badp is a list of the pdiff values that are less than peakdiffmin
    if (length(badp) != 0)
    {peaklist<-peaklist[-(badp+1)] # update peaklist, removing peaks identified by badp
    }
    
    # Do heart beat integration
    peakn=length(peaklist)
    rawdata$heartbeatcorrected_L <- 0
    rawdata$heartbeatcorrected_R <- 0
    for (p in 1:(peakn-1))
    {myrange=seq(peaklist[p],peaklist[p+1]) # the indices where the heartbeat will be replaced
    thisheart_L=mean(rawdata$normal_L[myrange]) # the new values that will be replaced
    thisheart_R=mean(rawdata$normal_R[myrange])
    rawdata$heartbeatcorrected_L[peaklist[p] : peaklist[p+1]]=thisheart_L
    rawdata$heartbeatcorrected_R[peaklist[p] : peaklist[p+1]]=thisheart_R
    if (p==1){
      rawdata$heartbeatcorrected_L[1:peaklist[p]] <- thisheart_L
      rawdata$heartbeatcorrected_R[1:peaklist[p]] <- thisheart_R
    }
    if (p==peakn-1){
      rawdata$heartbeatcorrected_L[peaklist[p] : mylen] <- thisheart_L
      rawdata$heartbeatcorrected_R[peaklist[p] : mylen] <- thisheart_R
    }
    }
    
    #----------------------------------------------------------
    # Save processed file in csv format
    mynewfile <- paste0(strsplit(myfile, '*.exp'), '_processed.csv')
    write.csv(rawdata, mynewfile, row.names=F)
    
    #----------------------------------------------------------
    ts_summary_data[j,] <- c(strsplit(myfile, '*.exp'), max(rawdata2$heartbeatcorrected_L), mean(rawdata2$heartbeatcorrected_L), max(rawdata2$heartbeatcorrected_R), mean(rawdata2$heartbeatcorrected_R))
    
    #---------------------------------------------------------------------------------------------------------------#
    
     myseq<-seq(1,length(rawdata[,1]),by=25)
     rawdata2<-rawdata[myseq,]
     
     blockends1<-cumsum(rle(rawdata2$stim1_on)$lengths)
     blockstarts1<-c(1,(blockends1+1)[-length(blockends1)])
     
     blockends2<-cumsum(rle(rawdata2$stim2_on)$lengths)
     blockstarts2<-c(1,(blockends2+1)[-length(blockends2)])
     
     rawdata2$epoch1<-rawdata2$epoch2<-rawdata2$time1<-rawdata2$time2<-rawdata2$adj_L1<-rawdata2$adj_R1<-rawdata2$adj_L2<-rawdata2$adj_R2<-rep(NA,length(rawdata2[,1]))

    for(i in 1:length(blockends))
    {
      rawdata2$epoch1[blockstarts1[i]:blockends1[i]]<-i
      rawdata2$epoch2[blockstarts2[i]:blockends2[i]]<-i

      rawdata2$time1[blockstarts1[i]:blockends1[i]]<-1:length(blockstarts1[i]:blockends1[i])
      rawdata2$time2[blockstarts2[i]:blockends2[i]]<-1:length(blockstarts2[i]:blockends2[i])

      rawdata2$adj_L1[blockstarts1[i]:blockends1[i]]<-rawdata2$heartbeatcorrected_L[blockstarts1[i]:blockends1[i]]-rawdata2$heartbeatcorrected_L[blockstarts1[i]]
      rawdata2$adj_L2[blockstarts2[i]:blockends2[i]]<-rawdata2$heartbeatcorrected_L[blockstarts2[i]:blockends2[i]]-rawdata2$heartbeatcorrected_L[blockstarts2[i]]

      rawdata2$adj_R1[blockstarts1[i]:blockends1[i]]<-rawdata2$heartbeatcorrected_R[blockstarts1[i]:blockends1[i]]-rawdata2$heartbeatcorrected_R[blockstarts1[i]]
      rawdata2$adj_R2[blockstarts2[i]:blockends2[i]]<-rawdata2$heartbeatcorrected_R[blockstarts2[i]:blockends2[i]]-rawdata2$heartbeatcorrected_R[blockstarts2[i]]

      #rawdata2$adj_L_c[blockstarts[i]:blockends[i]]<-rawdata2$heartbeatcorrected_L[blockstarts[i]:blockends[i]]-mean(rawdata2$heartbeatcorrected_L[blockstarts[i]:blockends[i]])
    }

     rawdata3a<-rawdata2 %>% filter(stim1_on == 1)
     rawdata3b<-rawdata2 %>% filter(stim2_on == 1)
     rawdata3a$epoch1<-factor(rawdata3a$epoch1)
     rawdata3a$epoch2<-factor(rawdata3a$epoch2)
     rawdata3b$epoch1<-factor(rawdata3b$epoch1)
     rawdata3b$epoch2<-factor(rawdata3b$epoch2)
     
    # 
    # #pdf(file = paste0(strsplit(myfile,'[.]')[[1]][1],'_HRF_plot_LEFT.pdf'))
     g1a<-ggplot(rawdata3a,aes(y=adj_L1,x=time1,colour=epoch1))+geom_line(show.legend = FALSE)+theme_bw()+stat_summary(fun.y = mean,geom = "line",colour="black")+stat_summary(fun.data = mean_cl_boot,geom = "ribbon",colour='grey',alpha=0.2)
     
     g1b<-ggplot(rawdata3b,aes(y=adj_L2,x=time2,colour=epoch2))+geom_line(show.legend = FALSE)+theme_bw()+stat_summary(fun.y = mean,geom = "line",colour="black")+stat_summary(fun.data = mean_cl_boot,geom = "ribbon",colour='grey',alpha=0.2)
    # #dev.off()
    # 
    # #pdf(file = paste0(strsplit(myfile,'[.]')[[1]][1],'_HRF_plot_RIGHT.pdf'))
     g2a<-ggplot(rawdata3a,aes(y=adj_R1,x=time1,colour=epoch1))+geom_line(show.legend = FALSE)+theme_bw()+stat_summary(fun.y = mean,geom = "line",colour="black")+stat_summary(fun.data = mean_cl_boot,geom = "ribbon",colour='grey',alpha=0.2)
     
     g2b<-ggplot(rawdata3b,aes(y=adj_R2,x=time2,colour=epoch2))+geom_line(show.legend = FALSE)+theme_bw()+stat_summary(fun.y = mean,geom = "line",colour="black")+stat_summary(fun.data = mean_cl_boot,geom = "ribbon",colour='grey',alpha=0.2)
    # #dev.off()
    
    
    
    
    #---------------------------------------------------------------------------------------------------------------#
    

    gamma1 = fmri.stimulus(scans = dim(rawdata)[1], onsets = c(1,1+which(diff(rawdata$stim1_on)!=0))[seq(2, length(c(1,1+which(diff(rawdata$stim1_on)!=0))), by = 2)], durations = 375, TR = 1/25,type="gamma",scale=1)
    
    gamma2 = fmri.stimulus(scans = dim(rawdata)[1], onsets = c(1,1+which(diff(rawdata$stim2_on)!=0))[seq(2, length(c(1,1+which(diff(rawdata$stim2_on)!=0))), by = 2)], durations = 125, TR = 1/25,type="gamma",scale=1)
    
    gamma = as.matrix(cbind(gamma1,gamma2))
    
    my_des<-fmri.design(gamma, order = order)
    
    #Left signal glm
    myfitL<-glm.fit(x=my_des,y=rawdata$heartbeatcorrected_L[c(seq(from=0, to=length(rawdata[,1]), by=25))[-1]],family=gaussian())
    
    class(myfitL) <- c(myfitL$class, c("glm", "lm"))
    
    #Right signal glm
    myfitR<-glm.fit(x=my_des,y=rawdata$heartbeatcorrected_R[c(seq(from=0, to=length(rawdata[,1]), by=25))[-1]],family=gaussian())
    
    class(myfitR) <- c(myfitR$class, c("glm", "lm"))
    
    
    glm.data[j,1] <- strsplit(basename(myfile),'[.]')[[1]][1]
    
    glm.data[j,(((order+3)*1)+2)] <- "gamma"
    glm.data[j,(((order+3)*2)+3)] <- "gamma"
    
    glm.data[j,2:(((order+3)*1)+1)] <- myfitL$coefficients
    glm.data[j,(((order+3)*1)+3):(((order+3)*2)+2)] <- myfitR$coefficients
    
    #Left
    myfit.diagL<-glm.diag(myfitL)
    # pdf(file = paste0(strsplit(myfile,'[.]')[[1]][1],'_diag_plot_LEFT.pdf'))
    # glm.diag.plots(switch(mychoiceL,boxcar=myfit1L,canonical=myfit2L,gamma=myfit3L), myfit.diagL)
    # dev.off()
    
    #Right
    myfit.diagR<-glm.diag(myfitR)
    # pdf(file = paste0(strsplit(myfile,'[.]')[[1]][1],'_diag_plot_RIGHT.pdf'))
    # glm.diag.plots(switch(mychoiceR,boxcar=myfit1R,canonical=myfit2R,gamma=myfit3R), myfit.diagR)
    # dev.off()
    
    myplotdat<-data.frame(y=c(rawdata$heartbeatcorrected_L[c(seq(from=0, to=length(rawdata[,1]), by=25))[-1]],rawdata$heartbeatcorrected_R[c(seq(from=0, to=length(rawdata[,1]), by=25))[-1]]),
                          x=c(rawdata$sec[c(seq(from=0, to=length(rawdata[,1]), by=25))[-1]],rawdata$sec[c(seq(from=0, to=length(rawdata[,1]), by=25))[-1]]),
                          fitted=c(myfitL$fitted.values,myfitR$fitted.values),Signal=rep(c("Left","Right"),each=length(rawdata$sec[c(seq(from=0, to=length(rawdata[,1]), by=25))[-1]])))


    g3<-ggplot(myplotdat,aes(y=y,x=x))+geom_point(colour='grey',alpha=0.5)+geom_line(aes(y=fitted),color="red")+theme_bw() + facet_grid(~Signal)

    # 
    # 
     g4<-ggarrange(g3,                                                 # First row with scatter plot
                ggarrange(g1a, g1b, g2a, g2b, ncol = 2,nrow=2, labels = c("B", "C","D", "E")), # Second row with box and dot plots
                   nrow = 2,
                   labels = "A"                                        # Labels of the scatter plot
     )

    #g4<-ggarrange(g1a, g1b, g2a, g2b, ncol = 2,nrow=2, labels = c("B", "C","D", "E"))

    g4<-annotate_figure(g4, top = text_grob(strsplit(myfile,'[.]')[[1]][1], face = "bold", size = 14))

    #pdf(file = paste0(strsplit(myfile,'[.]')[[1]][1],'_HRF_plot_signals.pdf'))

    #pdf(file = 'HRF_signals_plots.pdf', onefile = TRUE)
    #print(g1a)
    #print(g2a)
    #print(g1b)
    #print(g2b)
    #print(g3)
    print(g4)
    #dev.off()

    
    glm_data<-glm.data
  }
  write.csv(ts_summary_data,"ts_summary_data.csv")
  return(glm_data)    
}

#-----------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------#
#Set the order
order=3 #polynomial drift terms (2=quadratic, 3=cubic, etc...)
pdf(file = 'HRF_signals_plots.pdf', onefile = TRUE)
my_results<-fTCD_glm3(path=getwd(),order=order)
dev.off()
my_results_diff<-data.frame(matrix(NA,nrow=length(my_results[,1]),ncol=((order+3)+1)))
names(my_results_diff)<-c('ID',paste0('Dparam',(1:(order+3))))
my_results_diff$ID<-my_results$ID
for(i in 1:(order+3))
{
  my_results_diff[,(i+1)] <- my_results[,(i+1)]-my_results[,((order+3)+2+i)]
}
#-----------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------#


mylong_results_L<-gather(my_results[,1:(order+3+2)],key='param',value='beta',-c(ID,L_HRF))
mylong_results_R<-gather(my_results[,c(1,(((order+3)*1)+3):(((order+3)*2)+3))],key='param',value='beta',-c(ID,R_HRF))

mylong_results_D<-gather(my_results_diff[,1:(order+3+1)],key='param',value='beta',-ID)

names(mylong_results_L)[2]<-"HRF"
names(mylong_results_R)[2]<-"HRF"

mylong_results_L$Signal<-rep('Left',length(mylong_results_L$param))
mylong_results_R$Signal<-rep('Right',length(mylong_results_R$param))

#-----------------------------------------------------------------------------------------------------------------------#

mylong_results<-dplyr::bind_rows(mylong_results_L,mylong_results_R)

mylong_results$param<-substring(mylong_results$param, 2)

ggplot(mylong_results,aes(x=beta,group=Signal))+geom_density(aes(fill=Signal),alpha=0.5)+facet_wrap(~param,scales='free')+theme_bw()

ggplot(mylong_results,aes(x=HRF,group=Signal))+geom_bar(aes(fill=Signal,colour=Signal),alpha=0.5,position="dodge")+theme_bw()

#-----------------------------------------------------------------------------------------------------------------------#

ggplot(mylong_results_D,aes(x=beta))+geom_density(fill='blue',alpha=0.5)+facet_wrap(~param,scales='free')+theme_bw()

#-----------------------------------------------------------------------------------------------------------------------#


library(psych)
library(tidyverse)

#load LI based on old doppler analysis method

old_res<-read.csv("WordGen_results.csv")

old_res<-old_res %>% rename(ID=Filename)

compare_results<-merge(my_results_diff,old_res,by='ID',all.x = T)

#compare_results$New_LI <- ifelse(compare_results$Dparam1>0,)


fmri_data <- read.csv('/Volumes/PSYHOME/PSYRES/pthompson/DVMB/bruckett_reanalysis/Chapter5_fMRI_data.csv')

# Identify factors
factor_variables <- c('group_cat', 'group_lat', 'sex', 'hand_self_report', 'hand_QHP_cat', 'hand_EHI_cat', 
                      'Old_fTCD_wg_cat', 'Old_fTCD_pptt_cat', 'fTCD_wg_cat', 'fTCD_pptt_cat')
for (i in 1:length(factor_variables))
{factor_ind <- str_which(colnames(fmri_data), paste0('^',factor_variables[i]))
fmri_data[,factor_ind] <- as.factor(fmri_data[,factor_ind])}

# Relabel group_cat and sex factors for clarity
# NB: group_cat 0=typical; 1=atypical
# sex 0=male; 1=female
levels(fmri_data$group_cat) <- c('T', 'A')
levels(fmri_data$sex) <- c('M', 'F')


fmri_data<-fmri_data[,c('ID','fMRI_diff_wg_frontal','fMRI_diff_wg_temporal','fMRI_diff_wg_MCA')]

compare_results2<-merge(compare_results,fmri_data,by='ID')

psych::pairs.panels(compare_results2[,c('fMRI_diff_wg_frontal','fMRI_diff_wg_temporal','fMRI_diff_wg_MCA','LI','Dparam1')])
