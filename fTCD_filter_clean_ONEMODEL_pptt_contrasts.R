
# install.packages("remotes")
#remotes::install_github("centerforopenscience/osfr")


## Packages
library(osfr)
library(lme4)
library(utils)
require(dplyr)
require(tidyverse)
require(fmri)
require(boot)
require(ggpubr)
library(psych)

#---------------------------------------------------------------------------------------------------------------#

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
    osf_download() %>% unzip('Chpt4_fTCD_PPTT_rawdata.zip')
}
#---------------------------------------------------------------------------------------------------------------#

#unzip('Chpt4_fTCD_PPTT_rawdata.zip')

fTCD_glm4<-function(path,order)
{
  
  
  #---------------------------------------------------------------------------------------------------------------#
  
  filename1<-list.files(path,pattern = '.exp')
  ## Set parameters
  samplingrate <- 25 # Sampling rate after downsampling. Raw data is 100Hz, we take 1 in every 4 samples
  heartratemax <- 125
  
  
  glm.data<-data.frame(matrix(NA,nrow=length(filename1),ncol=(((order+4))+2)))
  names(glm.data)<-c('ID',paste0('param',(1:(order+4))),'HRF')
  

  for(j in 1:length(filename1))
  {
    print(filename1[j])
    
    #########################################
    # fTCD preprocessing for GLM analysis   #
    #                                       #
    # Created by z.woodhead 30th July 2019  #
    #########################################
    
    # This script takes a raw .exp datafile and preprocesses it ready for GLM analysis:
    #   - It creates a box car function showing when the task was ON or OFF
    #   - It normalises the fTCD signal to a mean of 100
    #   - It performs heart beat integration
    #   - It saves the processed data into a .csv file
    
    ## Read in raw data
    
    myfile <- filename1[j]
    
    mydata<-read.table(paste0(path,"/",myfile), skip = 6,  header =FALSE, sep ='\t')
    
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
    
    # Stimulus timings: PPTT starts 5 seconds after marker and continues for 20 seconds 
    # This includes 8 individual events, each 2.5 seconds long

    stim1_delay_sec <- 5
    stim1_delay_samples <- stim1_delay_sec * samplingrate
    stim1_length_sec <- 20
    stim1_length_samples <- stim1_length_sec * samplingrate
    
    # OPTION! If we want a more complicated model we can specify every event (i.e. all 8 stimuli within each block)
    
    # stim1_delay_sec <- 5
    # stim1_delay_samples <- stim1_delay_sec * samplingrate
    # stim1_event_length_sec <- .1    # to model all 8 events, we will specify the onset of each event only, with a 100ms boxcar
    # stim1_event_length_samples <- stim1_event_length_sec * samplingrate
    # stim1_event_interval_sec <- 2.5
    # stim1_n_events <- 8
    
    # stim2 not required for pptt
    
    # stim2_delay_sec <- 20
    # stim2_delay_samples <- stim2_delay_sec * samplingrate
    # stim2_length_sec <- 5
    # stim2_length_samples <- stim2_length_sec * samplingrate
    
    rest_length_sec <- 30
    rest_length_samples <- rest_length_sec * samplingrate
    
    rawdata$stim1_on <- 0
    # rawdata$stim2_on <- 0
    for (m in 1:norigmarkers){
      rawdata$stim1_on[(origmarkerlist[m]+stim1_delay_samples):(origmarkerlist[m]+stim1_delay_samples+stim1_length_samples)] <- 1
      
      # OPTION: here's the more complicated model:
      # 
      # for (i in 1:stim1_n_events){
      #   start_ind <- origmarkerlist[m] + stim1_delay_samples + stim1_event_interval_sec*(i-1)*samplingrate
      #   rawdata$stim1_on[start_ind : (start_ind + stim1_event_length_samples)] <- 1
      #   }

      # stim2 is not in use for pptt
      # rawdata$stim2_on[(origmarkerlist[m]+stim2_delay_samples):(origmarkerlist[m]+stim2_delay_samples+stim2_length_samples)] <- 1
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
    
    #---------------------------------------------------------------------------------------------------------------#
    # Save processed file in csv format
    mynewfile <- paste0(getwd(),"/Chpt4_fTCD_PPTT_rawdata/",strsplit(myfile, '*.exp'), '_processed.csv')
    write.csv(rawdata, mynewfile, row.names=F)

    #---------------------------------------------------------------------------------------------------------------#
    
    myseq<-seq(1,length(rawdata[,1]),by=25)
    rawdata2<-rawdata[myseq,]
    
    
    blockends1<-cumsum(rle(rawdata2$stim1_on)$lengths)
    blockstarts1<-c(1,(blockends1+1)[-length(blockends1)])
    
    # blockends2<-cumsum(rle(rawdata2$stim2_on)$lengths)
    # blockstarts2<-c(1,(blockends2+1)[-length(blockends2)])
    
    rawdata2$epoch1<-rawdata2$epoch2<-rawdata2$time1<-rawdata2$time2<-rawdata2$adj_L1<-rawdata2$adj_R1<-rawdata2$adj_L2<-rawdata2$adj_R2<-rep(NA,length(rawdata2[,1]))
    
    for(i in 1:length(blockends1))
    {
      rawdata2$epoch1[blockstarts1[i]:blockends1[i]]<-i
      # rawdata2$epoch2[blockstarts2[i]:blockends2[i]]<-i
      
      rawdata2$time1[blockstarts1[i]:blockends1[i]]<-1:length(blockstarts1[i]:blockends1[i])
      # rawdata2$time2[blockstarts2[i]:blockends2[i]]<-1:length(blockstarts2[i]:blockends2[i])
      
      rawdata2$adj_L1[blockstarts1[i]:blockends1[i]]<-rawdata2$heartbeatcorrected_L[blockstarts1[i]:blockends1[i]]-rawdata2$heartbeatcorrected_L[blockstarts1[i]]
      # rawdata2$adj_L2[blockstarts2[i]:blockends2[i]]<-rawdata2$heartbeatcorrected_L[blockstarts2[i]:blockends2[i]]-rawdata2$heartbeatcorrected_L[blockstarts2[i]]
      
      rawdata2$adj_R1[blockstarts1[i]:blockends1[i]]<-rawdata2$heartbeatcorrected_R[blockstarts1[i]:blockends1[i]]-rawdata2$heartbeatcorrected_R[blockstarts1[i]]
      # rawdata2$adj_R2[blockstarts2[i]:blockends2[i]]<-rawdata2$heartbeatcorrected_R[blockstarts2[i]:blockends2[i]]-rawdata2$heartbeatcorrected_R[blockstarts2[i]]
      
    }
    
    rawdata3a<-rawdata2 %>% filter(stim1_on == 1)
    # rawdata3b<-rawdata2 %>% filter(stim2_on == 1)
    rawdata3a$epoch1<-factor(rawdata3a$epoch1)
    rawdata3a$epoch2<-factor(rawdata3a$epoch2)
    # rawdata3b$epoch1<-factor(rawdata3b$epoch1)
    # rawdata3b$epoch2<-factor(rawdata3b$epoch2)
    
    #---------------------------------------------------------------------------------------------------------------#
  
    g1a<-ggplot(rawdata3a,aes(y=adj_L1,x=time1,colour=epoch1))+geom_line(show.legend = FALSE)+theme_bw()+stat_summary(fun.y = mean,geom = "line",colour="black")+stat_summary(fun.data = mean_cl_boot,geom = "ribbon",colour='grey',alpha=0.2)
    
    # g1b<-ggplot(rawdata3b,aes(y=adj_L2,x=time2,colour=epoch2))+geom_line(show.legend = FALSE)+theme_bw()+stat_summary(fun.y = mean,geom = "line",colour="black")+stat_summary(fun.data = mean_cl_boot,geom = "ribbon",colour='grey',alpha=0.2)

    g2a<-ggplot(rawdata3a,aes(y=adj_R1,x=time1,colour=epoch1))+geom_line(show.legend = FALSE)+theme_bw()+stat_summary(fun.y = mean,geom = "line",colour="black")+stat_summary(fun.data = mean_cl_boot,geom = "ribbon",colour='grey',alpha=0.2)
    
    # g2b<-ggplot(rawdata3b,aes(y=adj_R2,x=time2,colour=epoch2))+geom_line(show.legend = FALSE)+theme_bw()+stat_summary(fun.y = mean,geom = "line",colour="black")+stat_summary(fun.data = mean_cl_boot,geom = "ribbon",colour='grey',alpha=0.2)
  
    #---------------------------------------------------------------------------------------------------------------#
    
    fmri.stimulus.PT2<- function(scans = dim(rawdata)[1], onsets = c(1,1+which(diff(rawdata$stim1_on)!=0))[seq(2, length(c(1,1+which(diff(rawdata$stim1_on)!=0))), by = 2)], durations = 375, TR = 1/25,scale=1)
    {
      
      onsets <- onsets * TR
      durations <- durations * TR
      onsets <- onsets * scale
      durations <- durations * scale
      scans <- scans * TR * scale
      TR <- TR/scale
      no <- length(onsets)
      durations <- rep(durations, no)
      
      stimulus <- rep(0, ceiling(scans))
      for (i in 1:no) stimulus[onsets[i]:(onsets[i] + durations[i] - 
                                            1)] <- 1
      
      .gammaHRF <- function(t, par = NULL) {
        th <- 0.242 * par[1]
        1/(th * factorial(3)) * (t/th)^3 * exp(-t/th)
      }
      
      
      par <- 4
      
      y <- .gammaHRF(0:(28 * scale)/scale, par)
      
      stimulus <- convolve(stimulus, rev(y), type = "open")
      stimulus <- stimulus[unique((scale:scans)%/%(scale^2 * TR)) * scale^2 * TR]/(scale^2 * TR)
      stimulus <- stimulus - mean(stimulus)
      return(stimulus)
    }  
    
    #---------------------------------------------------------------------------------------------------------------#
    
    fmri.stimulus.PT3 <- function (scans = dim(rawdata)[1], onsets = c(1,1+which(diff(rawdata$stim1_on)!=0))[seq(2, length(c(1,1+which(diff(rawdata$stim1_on)!=0))), by = 2)], durations = 375, TR = 1/25,scale=1) 
    {
      
      onsets <- onsets * TR
      durations <- durations * TR
      
      onsets <- onsets * scale
      durations <- durations * scale
      scans <- scans * TR * scale
      TR <- TR/scale
      
      no <- length(onsets)
      
      durations <- rep(durations, no)
      
      stimulus <- rep(0, ceiling(scans))
      for (i in 1:no) {stimulus[onsets[i]:(onsets[i] + durations[i] - 1)] <- 1}
      
      y <- scale
      
      stimulus <- convolve(stimulus, rev(y), type = "open")
      stimulus <- stimulus[unique((scale:scans)%/%(scale^2 * TR)) * scale^2 * TR]/(scale^2 * TR)
      stimulus <- stimulus - mean(stimulus)
      
      return(stimulus)
    }
    
    #---------------------------------------------------------------------------------------------------------------#
    
    gamma1 = fmri.stimulus.PT2(scans = dim(rawdata)[1], onsets = c(1,1+which(diff(rawdata$stim1_on)!=0))[seq(2, length(c(1,1+which(diff(rawdata$stim1_on)!=0))), by = 2)], durations = 500, TR = 1/25,scale=1)

    gamma = c(gamma1,gamma1)
    
    my_des<-fmri.design(gamma, order = order)

    my_des<-cbind(my_des,rep(1:0,each=length(gamma1)))
    
    #---------------------------------------------------------------------------------------------------------------#
    
    mydata<-data.frame(y=c(rawdata$heartbeatcorrected_L[c(seq(from=0, to=length(rawdata[,1]), by=25))[-1]],rawdata$heartbeatcorrected_R[c(seq(from=0, to=length(rawdata[,1]), by=25))[-1]]),stim1=my_des[,1],t=my_des[,3],signal=as.factor(my_des[,6]),interaction=my_des[,1]*my_des[,6])
    
   
    myfit <- glm(y~stim1+t+I(t^2)+I(t^3)+signal+interaction,data=mydata,contrasts=c(0,1,0,0,0,0,0))

    names(myfit$coefficients)<-c("intercept","stim1","t","t_sqr","t_cub","signal","interaction")

    #---------------------------------------------------------------------------------------------------------------#
    
    glm.data[j,1] <- strsplit(basename(myfile),'[.]')[[1]][1]
    
    glm.data[j,(((order+4)*1)+2)] <- "gamma"
    
    glm.data[j,2:(((order+4)*1)+1)] <- myfit$coefficients
    
    #---------------------------------------------------------------------------------------------------------------#
    
    pframe<-with(rawdata,expand.grid(t=seq(min(sec),max(sec),length=length(rawdata$heartbeatcorrected_L[c(seq(from=0, to=length(rawdata[,1]), by=25))[-1]])),signal=c(0,1)))
    
    pframe<-data.frame(stim1=c(gamma1,gamma1),t=pframe[,1],t_sqr=(pframe[,1])^2,t_cub=(pframe[,1])^3,signal=pframe[,2])
    
    myplotdat<-data.frame(y=c(rawdata$heartbeatcorrected_L[c(seq(from=0, to=length(rawdata[,1]), by=25))[-1]],rawdata$heartbeatcorrected_R[c(seq(from=0, to=length(rawdata[,1]), by=25))[-1]]),
                          x=c(rawdata$sec[c(seq(from=0, to=length(rawdata[,1]), by=25))[-1]],rawdata$sec[c(seq(from=0, to=length(rawdata[,1]), by=25))[-1]]),
                          fitted=predict(myfit),Signal=rep(c("Left","Right"),each=length(rawdata$sec[c(seq(from=0, to=length(rawdata[,1]), by=25))[-1]])))
    
    
    
    
    g3<-ggplot(myplotdat,aes(y=y,x=x,colour=Signal))+geom_point(colour='grey',alpha=0.5)+geom_line(aes(y=fitted))+theme_bw()
    
    g4<-ggarrange(g3, ggarrange(g1a, g1b, g2a, g2b, ncol = 2,nrow=2, labels = c("B", "C","D", "E")), nrow = 2,labels = "A")

    g4<-annotate_figure(g4, top = text_grob(strsplit(myfile,'[.]')[[1]][1], face = "bold", size = 14))
    
    print(g4)
    
    glm_data<-glm.data
  }

  return(glm_data)    
}

#-----------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------#
#Set the order
order=3 #polynomial drift terms (2=quadratic, 3=cubic, etc...)
pdf(file = 'HRF_signals_plots_PPTT.pdf', onefile = TRUE)
my_results<-fTCD_glm4(path=paste0(getwd(),"/Chpt4_fTCD_PPTT_rawdata"),order=order)
dev.off()
#stop()
#-----------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------#


mylong_results<-gather(my_results,key='param',value='beta',-c(ID,HRF))

names(mylong_results)[2]<-"HRF"


ggplot(mylong_results,aes(x=beta))+geom_density(fill='blue',alpha=0.5)+facet_wrap(~param,scales='free')+theme_bw()



#-----------------------------------------------------------------------------------------------------------------------#

#load LI based on old doppler analysis method

old_res<-read.csv(paste0(getwd(),"/Chpt4_fTCD_PPTT_rawdata/","PPTT_results.csv"))

old_res<-old_res %>% rename(ID=Filename)

compare_results<-merge(my_results,old_res,by='ID',all.x = T)

#compare_results$New_LI <- ifelse(compare_results$Dparam1>0,)


fmri_data <- read.csv('/Volumes/PSYHOME/PSYRES/pthompson/DVMB/bruckett_reanalysis/Chapter5_fMRI_data.csv')


fmri_data<-fmri_data[,c('ID','fMRI_diff_pptt_frontal','fMRI_diff_pptt_temporal','fMRI_diff_pptt_MCA')]

compare_results$ID<-substring(compare_results$ID,1,6)
fmri_data$ID<-substring(fmri_data$ID,1,6)

compare_results2<-merge(compare_results,fmri_data,by='ID')

psych::pairs.panels(compare_results2[,c('fMRI_diff_pptt_frontal','fMRI_diff_pptt_temporal','fMRI_diff_pptt_MCA','LI','param1','param2','param3','param4','param5','param6','param7')])


