
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
  
  osf_retrieve_file("https://osf.io/j62x4/") %>%
    osf_download() %>% unzip('Holly_fTCD_data_run1.zip')
}
#---------------------------------------------------------------------------------------------------------------#

#unzip('Holly_fTCD_data_run1.zip')

fTCD_glm_multi<-function(path,order)
{
  
  
  #---------------------------------------------------------------------------------------------------------------#
  
  filename1<-list.files(path,pattern = '.exp')
  ## Set parameters
  samplingrate <- 25 # Sampling rate after downsampling. Raw data is 100Hz, we take 1 in every 4 samples
  heartratemax <- 125
  
  
  glm.data<-data.frame(matrix(NA,nrow=length(filename1),ncol=(((order+9))+2)))
  names(glm.data)<-c('ID',paste0('param',(1:(order+9))),'HRF')
  
  
  for(j in 1:length(filename1))
  {
    print(filename1[j])
    
    #########################################
    # fTCD preprocessing for GLM analysis   #
    #                                       #
    # Created by z.woodhead 30th July 2019  #
    # Edited  by z. woodhead 3rd Oct 2019   #
    #########################################
    
    # This script takes a raw .exp datafile and preprocesses it ready for GLM analysis:
    #   - It creates a box car function showing when the task was ON or OFF
    #   - It normalises the fTCD signal to a mean of 100
    #   - It performs heart beat integration
    #   - It saves the processed data into a .csv file
    
    ## Read in raw data
    
    myfile <- filename1[j]
    
    mydata<-read.table(paste0(path,"/",myfile), skip = 6,  header =FALSE, sep ='\t')
    
    wantcols = c(2,3,4,9) #sec, L, R,marker #select columns of interest to put in shortdat
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
    
    # Stimulus order: In this task, there were three tasks:
        # 1) Word Generation (WG)
        # 2) Sentence Generation (SG)
        # 3) List Generation (LG)
    
    # Here's the list of the order they appeared in run 1:
    stim_order <- c(1,2,3,2,3,1,3,1,2,3,2,1,2,1,3,1,3,2,1,2,3,2,3,1,3,1,2,3,2,1)
    
    # Stimulus timings: Each task has the same timings. Each task is comprised of two stimuli: stim1 = covert speech generation; stim2 = overt speech generation
    # The stimulus starts 3 seconds after marker (stim1_delay_sec = 3)
    # Covert generation lasts for 12 seconds (stim1_length_sec = 12)
    # Overt generation ('reporting') starts immediately after, i.e. 15 seconds after marker (stim2_delay_sec = 15
    # Overt generation lasts for 5 seconds (stim2_length_sec = 5)

    stim1_delay_sec <- 3
    stim1_delay_samples <- stim1_delay_sec * samplingrate
    stim1_length_sec <- 12
    stim1_length_samples <- stim1_length_sec * samplingrate
    
    stim2_delay_sec <- 15
    stim2_delay_samples <- stim2_delay_sec * samplingrate
    stim2_length_sec <- 5
    stim2_length_samples <- stim2_length_sec * samplingrate 
    
    # There is 10 seconds of rest between trials
    rest_length_sec <- 10
    rest_length_samples <- rest_length_sec * samplingrate
    
    rawdata$WG_stim1_on <- 0
    rawdata$WG_stim2_on <- 0
    rawdata$SG_stim1_on <- 0
    rawdata$SG_stim2_on <- 0
    rawdata$LG_stim1_on <- 0
    rawdata$LG_stim2_on <- 0
    
    # rawdata$stim2_on <- 0
    for (m in 1:norigmarkers){
      mytask <- stim_order[m]
      
      if (mytask == 1) {
        rawdata$WG_stim1_on[(origmarkerlist[m]+stim1_delay_samples):(origmarkerlist[m]+stim1_delay_samples+stim1_length_samples)] <- 1
        rawdata$WG_stim2_on[(origmarkerlist[m]+stim2_delay_samples):(origmarkerlist[m]+stim2_delay_samples+stim2_length_samples)] <- 1
      }
     
      if (mytask == 2) {
        rawdata$SG_stim1_on[(origmarkerlist[m]+stim1_delay_samples):(origmarkerlist[m]+stim1_delay_samples+stim1_length_samples)] <- 1
        rawdata$SG_stim2_on[(origmarkerlist[m]+stim2_delay_samples):(origmarkerlist[m]+stim2_delay_samples+stim2_length_samples)] <- 1
      }
      
      if (mytask == 3) {
        rawdata$LG_stim1_on[(origmarkerlist[m]+stim1_delay_samples):(origmarkerlist[m]+stim1_delay_samples+stim1_length_samples)] <- 1
        rawdata$LG_stim2_on[(origmarkerlist[m]+stim2_delay_samples):(origmarkerlist[m]+stim2_delay_samples+stim2_length_samples)] <- 1
      }
      
    }
    
    # print(dim(rawdata))
    # print(
    # ts.plot(rawdata$WG_stim1_on, col='red')
    # lines(rawdata$WG_stim2_on, col='red')
    # lines(rawdata$SG_stim1_on, col='blue')
    # lines(rawdata$SG_stim2_on, col='blue')
    # lines(rawdata$LG_stim1_on, col='green')
    # lines(rawdata$LG_stim2_on, col='green')
    #)
    
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
      & (rawdata$L[i] > rawdata$L[i+4])
      & (rawdata$L[i] > rawdata$L[i+5]))
    {peaklist=c(peaklist,i)}
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
    mynewfile <- paste0(getwd(),"/Holly_fTCD_data_run1/",strsplit(myfile, '*.exp'), '_processed.csv')
    write.csv(rawdata, mynewfile, row.names=F)
    
    #---------------------------------------------------------------------------------------------------------------#
    
    #myseq<-seq(1,length(rawdata[,1]),by=25)
    rawdata2<-rawdata#[myseq,]
    #---------------------------------------------------------------------------------------------------------------#
    # WG
    WG_blockends1<-cumsum(rle(rawdata2$WG_stim1_on)$lengths)
    WG_blockstarts1<-c(1,(WG_blockends1+1)[-length(WG_blockends1)])
    
    WG_blockends2<-cumsum(rle(rawdata2$WG_stim2_on)$lengths)
    WG_blockstarts2<-c(1,(WG_blockends2+1)[-length(WG_blockends2)])
    
    rawdata2$WG_epoch1<-rawdata2$WG_epoch2<-rawdata2$WG_time1<-rawdata2$WG_time2<-rawdata2$WG_adj_L1<-rawdata2$WG_adj_R1<-rawdata2$WG_adj_L2<-rawdata2$WG_adj_R2<-rep(NA,length(rawdata2[,1]))
    
    for(i in 1:length(WG_blockends1))
    {
      rawdata2$WG_epoch1[WG_blockstarts1[i]:WG_blockends1[i]]<-i
      rawdata2$WG_epoch2[WG_blockstarts2[i]:WG_blockends2[i]]<-i
      
      rawdata2$WG_time1[WG_blockstarts1[i]:WG_blockends1[i]]<-1:length(WG_blockstarts1[i]:WG_blockends1[i])
      rawdata2$WG_time2[WG_blockstarts2[i]:WG_blockends2[i]]<-1:length(WG_blockstarts2[i]:WG_blockends2[i])
      
      rawdata2$WG_adj_L1[WG_blockstarts1[i]:WG_blockends1[i]]<-rawdata2$heartbeatcorrected_L[WG_blockstarts1[i]:WG_blockends1[i]]-rawdata2$heartbeatcorrected_L[WG_blockstarts1[i]]
      rawdata2$WG_adj_L2[WG_blockstarts2[i]:WG_blockends2[i]]<-rawdata2$heartbeatcorrected_L[WG_blockstarts2[i]:WG_blockends2[i]]-rawdata2$heartbeatcorrected_L[WG_blockstarts2[i]]
      
      rawdata2$WG_adj_R1[WG_blockstarts1[i]:WG_blockends1[i]]<-rawdata2$heartbeatcorrected_R[WG_blockstarts1[i]:WG_blockends1[i]]-rawdata2$heartbeatcorrected_R[WG_blockstarts1[i]]
      rawdata2$WG_adj_R2[WG_blockstarts2[i]:WG_blockends2[i]]<-rawdata2$heartbeatcorrected_R[WG_blockstarts2[i]:WG_blockends2[i]]-rawdata2$heartbeatcorrected_R[WG_blockstarts2[i]]
    }
    #---------------------------------------------------------------------------------------------------------------#
    #SG
    SG_blockends1<-cumsum(rle(rawdata2$SG_stim1_on)$lengths)
    SG_blockstarts1<-c(1,(SG_blockends1+1)[-length(SG_blockends1)])
    
    SG_blockends2<-cumsum(rle(rawdata2$SG_stim2_on)$lengths)
    SG_blockstarts2<-c(1,(SG_blockends2+1)[-length(SG_blockends2)])
    
    rawdata2$SG_epoch1<-rawdata2$SG_epoch2<-rawdata2$SG_time1<-rawdata2$SG_time2<-rawdata2$SG_adj_L1<-rawdata2$SG_adj_R1<-rawdata2$SG_adj_L2<-rawdata2$SG_adj_R2<-rep(NA,length(rawdata2[,1]))
    
    for(i in 1:length(SG_blockends1))
    {
      rawdata2$SG_epoch1[SG_blockstarts1[i]:SG_blockends1[i]]<-i
      rawdata2$SG_epoch2[SG_blockstarts2[i]:SG_blockends2[i]]<-i
      
      rawdata2$SG_time1[SG_blockstarts1[i]:SG_blockends1[i]]<-1:length(SG_blockstarts1[i]:SG_blockends1[i])
      rawdata2$SG_time2[SG_blockstarts2[i]:SG_blockends2[i]]<-1:length(SG_blockstarts2[i]:SG_blockends2[i])
      
      rawdata2$SG_adj_L1[SG_blockstarts1[i]:SG_blockends1[i]]<-rawdata2$heartbeatcorrected_L[SG_blockstarts1[i]:SG_blockends1[i]]-rawdata2$heartbeatcorrected_L[SG_blockstarts1[i]]
      rawdata2$SG_adj_L2[SG_blockstarts2[i]:SG_blockends2[i]]<-rawdata2$heartbeatcorrected_L[SG_blockstarts2[i]:SG_blockends2[i]]-rawdata2$heartbeatcorrected_L[SG_blockstarts2[i]]
      
      rawdata2$SG_adj_R1[SG_blockstarts1[i]:SG_blockends1[i]]<-rawdata2$heartbeatcorrected_R[SG_blockstarts1[i]:SG_blockends1[i]]-rawdata2$heartbeatcorrected_R[SG_blockstarts1[i]]
      rawdata2$SG_adj_R2[SG_blockstarts2[i]:SG_blockends2[i]]<-rawdata2$heartbeatcorrected_R[SG_blockstarts2[i]:SG_blockends2[i]]-rawdata2$heartbeatcorrected_R[SG_blockstarts2[i]]
    }
    #---------------------------------------------------------------------------------------------------------------#
    #LG
    LG_blockends1<-cumsum(rle(rawdata2$LG_stim1_on)$lengths)
    LG_blockstarts1<-c(1,(LG_blockends1+1)[-length(LG_blockends1)])
    
    LG_blockends2<-cumsum(rle(rawdata2$LG_stim2_on)$lengths)
    LG_blockstarts2<-c(1,(LG_blockends2+1)[-length(LG_blockends2)])
    
    rawdata2$LG_epoch1<-rawdata2$LG_epoch2<-rawdata2$LG_time1<-rawdata2$LG_time2<-rawdata2$LG_adj_L1<-rawdata2$LG_adj_R1<-rawdata2$LG_adj_L2<-rawdata2$LG_adj_R2<-rep(NA,length(rawdata2[,1]))
    
    for(i in 1:length(LG_blockends1))
    {
      rawdata2$LG_epoch1[LG_blockstarts1[i]:LG_blockends1[i]]<-i
      rawdata2$LG_epoch2[LG_blockstarts2[i]:LG_blockends2[i]]<-i
      
      rawdata2$LG_time1[LG_blockstarts1[i]:LG_blockends1[i]]<-1:length(LG_blockstarts1[i]:LG_blockends1[i])
      rawdata2$LG_time2[LG_blockstarts2[i]:LG_blockends2[i]]<-1:length(LG_blockstarts2[i]:LG_blockends2[i])
      
      rawdata2$LG_adj_L1[LG_blockstarts1[i]:LG_blockends1[i]]<-rawdata2$heartbeatcorrected_L[LG_blockstarts1[i]:LG_blockends1[i]]-rawdata2$heartbeatcorrected_L[LG_blockstarts1[i]]
      rawdata2$LG_adj_L2[LG_blockstarts2[i]:LG_blockends2[i]]<-rawdata2$heartbeatcorrected_L[LG_blockstarts2[i]:LG_blockends2[i]]-rawdata2$heartbeatcorrected_L[LG_blockstarts2[i]]
      
      rawdata2$LG_adj_R1[LG_blockstarts1[i]:LG_blockends1[i]]<-rawdata2$heartbeatcorrected_R[LG_blockstarts1[i]:LG_blockends1[i]]-rawdata2$heartbeatcorrected_R[LG_blockstarts1[i]]
      rawdata2$LG_adj_R2[LG_blockstarts2[i]:LG_blockends2[i]]<-rawdata2$heartbeatcorrected_R[LG_blockstarts2[i]:LG_blockends2[i]]-rawdata2$heartbeatcorrected_R[LG_blockstarts2[i]]
      
    }
    #---------------------------------------------------------------------------------------------------------------#
    rawdata3a<-rawdata2 %>% filter(stim1_on == 1)
    rawdata3b<-rawdata2 %>% filter(stim2_on == 1)
    rawdata3a$epoch1<-factor(rawdata3a$epoch1)
    rawdata3a$epoch2<-factor(rawdata3a$epoch2)
    rawdata3b$epoch1<-factor(rawdata3b$epoch1)
    rawdata3b$epoch2<-factor(rawdata3b$epoch2)
    
    #---------------------------------------------------------------------------------------------------------------#
    
    g1a<-ggplot(rawdata3a,aes(y=adj_L1,x=time1,colour=epoch1))+geom_line(show.legend = FALSE)+theme_bw()+stat_summary(fun.y = mean,geom = "line",colour="black")+stat_summary(fun.data = mean_cl_boot,geom = "ribbon",colour='grey',alpha=0.2)
    
    g2a<-ggplot(rawdata3a,aes(y=adj_R1,x=time1,colour=epoch1))+geom_line(show.legend = FALSE)+theme_bw()+stat_summary(fun.y = mean,geom = "line",colour="black")+stat_summary(fun.data = mean_cl_boot,geom = "ribbon",colour='grey',alpha=0.2)
    
    
    #---------------------------------------------------------------------------------------------------------------#
    
    fmri.stimulus.PT2<- function(scans = dim(rawdata)[1], onsets = c(1,1+which(diff(rawdata$stim1_on)!=0))[seq(2, length(c(1,1+which(diff(rawdata$stim1_on)!=0))), by = 2)], durations = stim1_event_length_samples, TR = 1/25,scale=1)
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
      
      
      par <- floor((durations[1]/28)*4)
      
      y <- .gammaHRF(0:(durations[1] * scale)/scale, par)
      
      stimulus <- convolve(stimulus, rev(y), type = "open")

      stimulus <- stimulus[unique((scale:scans)%/%(scale^2 * TR)) * scale^2 * TR]/(scale^2 * TR)
      stimulus <- stimulus - mean(stimulus)
      return(stimulus)
    }  
    

    #---------------------------------------------------------------------------------------------------------------#
    
    gamma1 = fmri.stimulus.PT2(scans = dim(rawdata)[1], onsets = c(1,1+which(diff(rawdata$WG_stim1_on)!=0)), durations = stim1_length_samples, TR = 1,scale=1)
    
    gamma2 = fmri.stimulus.PT2(scans = dim(rawdata)[1], onsets = c(1,1+which(diff(rawdata$WG_stim2_on)!=0)), durations = stim2_length_samples, TR = 1,scale=1)
    
    gamma3 = fmri.stimulus.PT2(scans = dim(rawdata)[1], onsets = c(1,1+which(diff(rawdata$SG_stim1_on)!=0)), durations = stim1_length_samples, TR = 1,scale=1)
    
    gamma4 = fmri.stimulus.PT2(scans = dim(rawdata)[1], onsets = c(1,1+which(diff(rawdata$SG_stim2_on)!=0)), durations = stim2_length_samples, TR = 1,scale=1)
    
    gamma5 = fmri.stimulus.PT2(scans = dim(rawdata)[1], onsets = c(1,1+which(diff(rawdata$LG_stim1_on)!=0)), durations = stim1_length_samples, TR = 1,scale=1)
    
    gamma6 = fmri.stimulus.PT2(scans = dim(rawdata)[1], onsets = c(1,1+which(diff(rawdata$LG_stim2_on)!=0)), durations = stim2_length_samples, TR = 1,scale=1)
    
    #---------------------------------------------------------------------------------------------------------------# 
    
    gamma = as.matrix(cbind(gamma1,gamma2,gamma3,gamma4,gamma5,gamma6))
    
    gamma = rbind(gamma,gamma)
    
    #---------------------------------------------------------------------------------------------------------------#
    
    my_des<-fmri.design(gamma, order = order)
    
    my_des<-cbind(my_des,rep(1:0,each=length(gamma1)))
    
    #---------------------------------------------------------------------------------------------------------------#
    
    mydata<-data.frame(y=c(rawdata$heartbeatcorrected_L,rawdata$heartbeatcorrected_R),stim1=my_des[,1],stim2=my_des[,2],stim3=my_des[,3],stim4=my_des[,4],stim5=my_des[,5],stim6=my_des[,6],t=my_des[,7],signal=as.factor(my_des[,11]),interaction=my_des[,1]*my_des[,11])
    
    
    myfit <- glm(y~stim1+stim2+stim3+stim4+stim5+stim6+t+I(t^2)+I(t^3)+signal+interaction,data=mydata,contrasts=c(0,0,1,0,-1,0,0,0,0,0,0,0))
    
    class(myfit) <- c(myfit$class, c("glm", "lm"))
    names(myfit$coefficients)<-c("stim1","stim2","stim3","stim4","stim5","stim6","intercept","t","t_sqr","t_cub","signal","interaction")
    
    #---------------------------------------------------------------------------------------------------------------#
    
    glm.data[j,1] <- strsplit(basename(myfile),'[.]')[[1]][1]
    
    glm.data[j,(((order+9)*1)+2)] <- "gamma"
    
    glm.data[j,2:(((order+9)*1)+1)] <- myfit$coefficients
    
    #---------------------------------------------------------------------------------------------------------------#
    
    pframe<-with(rawdata,expand.grid(t=seq(min(sec),max(sec),length=length(rawdata$heartbeatcorrected_L)),signal=c(0,1)))
    
    pframe<-data.frame(stim1=c(gamma1,gamma1),stim2=c(gamma2,gamma2),stim3=c(gamma3,gamma3),stim4=c(gamma4,gamma4),stim5=c(gamma5,gamma5),stim6=c(gamma6,gamma6),t=pframe[,1],t_sqr=(pframe[,1])^2,t_cub=(pframe[,1])^3,signal=pframe[,2])
    
    myplotdat<-data.frame(y=c(rawdata$heartbeatcorrected_L,rawdata$heartbeatcorrected_R),
                          x=c(rawdata$sec,rawdata$sec),
                          fitted=predict(myfit),Signal=rep(c("Left","Right"),each=length(rawdata$sec)))
    
    #---------------------------------------------------------------------------------------------------------------# 
    
    
    g3<-ggplot(myplotdat,aes(y=y,x=x,colour=Signal))+geom_point(colour='grey',alpha=0.5)+geom_line(aes(y=fitted))+theme_bw()
    
    g4<-ggarrange(g3, ggarrange(g1a, g2a, ncol = 2,nrow=1, labels = c("B", "C")), nrow = 2,labels = "A")
    
    g4<-annotate_figure(g4, top = text_grob(strsplit(myfile,'[.]')[[1]][1], face = "bold", size = 14))
    
    print(g4)
    
    glm_data<-glm.data
  }
  
  return(glm_data)    
}
#
#-----------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------#
#Set the order
order=3 #polynomial drift terms (2=quadratic, 3=cubic, etc...)
pdf(file = 'HRF_signals_plots_multi_stim_holly.pdf', onefile = TRUE)
my_results_multi<-fTCD_glm_multi(path=paste0(getwd(),"/Holly_fTCD_data_run1"),order=order)
dev.off()
#stop()
#-----------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------#


mylong_resultsM<-gather(my_results_multi,key='param',value='beta',-c(ID,HRF))

names(mylong_resultsM)[2]<-"HRF"


ggplot(mylong_resultsM,aes(x=beta))+geom_density(fill='blue',alpha=0.5)+facet_wrap(~param,scales='free')+theme_bw()



#-----------------------------------------------------------------------------------------------------------------------#
# 
# #load LI based on old doppler analysis method
# 
# old_res<-read.csv(paste0(getwd(),"/Chpt4_fTCD_PPTT_rawdata/","PPTT_results.csv"))
# 
# old_res<-old_res %>% rename(ID=Filename)
# 
# compare_resultsM<-merge(my_results_multi,old_res,by='ID',all.x = T)
# 
# 
# 
# 
# fmri_data <- read.csv('/Volumes/PSYHOME/PSYRES/pthompson/DVMB/bruckett_reanalysis/Chapter5_fMRI_data.csv')
# 
# 
# fmri_data<-fmri_data[,c('ID','fMRI_diff_pptt_frontal','fMRI_diff_pptt_temporal','fMRI_diff_pptt_MCA')]
# 
# compare_resultsM$ID<-substring(compare_results$ID,1,6)
# fmri_data$ID<-substring(fmri_data$ID,1,6)
# 
# compare_results2M<-merge(compare_resultsM,fmri_data,by='ID')
# 
# psych::pairs.panels(compare_results2M[,c('fMRI_diff_pptt_frontal','fMRI_diff_pptt_temporal','fMRI_diff_pptt_MCA','LI','param1','param2','param3','param4','param5','param6','param7','param8','param9','param10','param11','param12')])
# 
# 
