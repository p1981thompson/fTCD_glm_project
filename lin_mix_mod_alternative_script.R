
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
library(contrast)
library(lmerTest)
library(cladoRcpp) # turbo charges the convolve function used to convolve the stimulus and HRF (supper slow previously).
library(data.table)


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

fTCD_glm_multi_lmm<-function(path,order)
{
  
  
  #---------------------------------------------------------------------------------------------------------------#
  
  filename1<-list.files(path,pattern = '.exp')
  ## Set parameters
  samplingrate <- 25 # Sampling rate after downsampling. Raw data is 100Hz, we take 1 in every 4 samples
  heartratemax <- 125
  
  
  glm.data<-data.frame(matrix(NA,nrow=length(filename1),ncol=14))
  names(glm.data)<-c('ID',paste0('param',1:12),'p_value_contrast')
  
  Rawdata<-0
  
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
    rawdata$ID<-strsplit(myfile, '*.exp')
    
   Rawdata<-rbind(Rawdata,rawdata)
   
   processed_data<-Rawdata
  }    
    #---------------------------------------------------------------------------------------------------------------#
    # Save processed file in csv format
    # mynewfile <- paste0(getwd(),"/Holly_fTCD_data_run1/",strsplit(myfile, '*.exp'), '_processed.csv')
    # write.csv(rawdata, mynewfile, row.names=F)
    # 
    #---------------------------------------------------------------------------------------------------------------#
    
  
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
      
      stimulus <- rcpp_convolve(a=stimulus, b=rev(y))
      
      stimulus <- stimulus[unique((scale:scans)%/%(scale^2 * TR)) * scale^2 * TR]/(scale^2 * TR)
      stimulus <- stimulus - mean(stimulus)
      return(stimulus)
    }  
    
    #---------------------------------------------------------------------------------------------------------------#
    bigDes <- matrix(, nrow = 0, ncol = 11)

    for(j in unique(processed_data$ID))
    {
    gamma1 = fmri.stimulus.PT2(scans=dim(processed_data[processed_data$ID==j,])[1],onsets = c(1,1+which(diff(processed_data$WG_stim1_on)!=0)), durations = stim1_length_samples, TR = 1,scale=1)
 
    gamma2 = fmri.stimulus.PT2(scans=dim(processed_data[processed_data$ID==j,])[1], onsets = c(1,1+which(diff(processed_data$WG_stim2_on)!=0)), durations = stim2_length_samples, TR = 1,scale=1)

    gamma3 = fmri.stimulus.PT2(scans=dim(processed_data[processed_data$ID==j,])[1], onsets = c(1,1+which(diff(processed_data$SG_stim1_on)!=0)), durations = stim1_length_samples, TR = 1,scale=1)

    gamma4 = fmri.stimulus.PT2(scans=dim(processed_data[processed_data$ID==j,])[1], onsets = c(1,1+which(diff(processed_data$SG_stim2_on)!=0)), durations = stim2_length_samples, TR = 1,scale=1)

    gamma5 = fmri.stimulus.PT2(scans=dim(processed_data[processed_data$ID==j,])[1], onsets = c(1,1+which(diff(processed_data$LG_stim1_on)!=0)), durations = stim1_length_samples, TR = 1,scale=1)

    gamma6 = fmri.stimulus.PT2(scans=dim(processed_data[processed_data$ID==j,])[1], onsets = c(1,1+which(diff(processed_data$LG_stim2_on)!=0)), durations = stim2_length_samples, TR = 1,scale=1)

    #---------------------------------------------------------------------------------------------------------------#
    gamma = as.matrix(cbind(gamma1,gamma2,gamma3,gamma4,gamma5,gamma6))

    gamma = rbind(gamma,gamma)

    #---------------------------------------------------------------------------------------------------------------#
    
    my_des<-fmri.design(gamma, order = order)

    my_des<-cbind(my_des,rep(1:0,each=length(gamma1)))

    bigDes<-rbind(bigDes,my_des)

    }

    #---------------------------------------------------------------------------------------------------------------# 
    #print(class(bigDes))
    #print(length(c(processed_data$heartbeatcorrected_L,processed_data$heartbeatcorrected_R)))
    
    # added data.table rather than data.frame, but now can't evalulate without do call. it can't handle id variable, so need to look into that (9/10/19)
    
    mydata<-data.table(y=c(processed_data$heartbeatcorrected_L,processed_data$heartbeatcorrected_R),stim1=bigDes[,1],stim2=bigDes[,2],stim3=bigDes[,3],stim4=bigDes[,4],stim5=bigDes[,5],stim6=bigDes[,6],t=bigDes[,8],signal=as.factor(bigDes[,11]),stim3_signal=bigDes[,3]*bigDes[,11],stim5_signal=bigDes[,5]*bigDes[,11],id=processed_data$ID)
    #---------------------------------------------------------------------------------------------------------------#
    # contrasts reference: https://osf.io/6kudn/ and https://psyarxiv.com/crx4m/ 
    print(str(mydata))
    #contrasts
    mydata$stim3_signal_adj <- (mydata$stim3_signal+mydata$stim5_signal)
    mydata$stim5_signal_adj <- mydata$stim5_signal
    
    #---------------------------------------------------------------------------------------------------------------#
    
    myfit <-mydata %>% do(lmer(y~stim1+stim2+stim3+stim4+stim5+stim6+t+I(t^2)+I(t^3)+signal+stim3_signal_adj+stim5_signal_adj + (1+stim1+stim2+stim3+stim4+stim5+stim6+signal+stim3_signal_adj+stim5_signal_adj|id),.))
 
    #---------------------------------------------------------------------------------------------------------------#
    
    glm.data[j,1] <- strsplit(basename(myfile),'[.]')[[1]][1]
    
    glm.data[j,2:13] <- fixef(myfit)
    
    glm.data[j,14] <- summary(myfit)$coefficients[12,5]
    
    #---------------------------------------------------------------------------------------------------------------#
    
    pframe<-with(processed_data,expand.grid(t=seq(min(sec),max(sec),length=length(heartbeatcorrected_L)),signal=c(0,1)))
    
    pframe<-data.frame(stim1=c(gamma1,gamma1),stim2=c(gamma2,gamma2),stim3=c(gamma3,gamma3),stim4=c(gamma4,gamma4),stim5=c(gamma5,gamma5),stim6=c(gamma6,gamma6),t=pframe[,1],t_sqr=(pframe[,1])^2,t_cub=(pframe[,1])^3,signal=pframe[,2])
    
    myplotdat<-data.frame(y=c(processed_data$heartbeatcorrected_L,processed_data$heartbeatcorrected_R), x=processed_data$sec,
                          fitted=predict(myfit),Signal=rep(c("Left","Right"),each=length(rawdata$sec)))

    #---------------------------------------------------------------------------------------------------------------# 
    
    
    g3<-ggplot(myplotdat,aes(y=y,x=x,colour=Signal))+geom_point(colour='grey',alpha=0.5)+geom_line(aes(y=fitted))+theme_bw()
    
    print(g3)
    
    glm_data<-glm.data

  
  return(glm_data)    
}

#
#-----------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------#
#Set the order
order=3 #polynomial drift terms (2=quadratic, 3=cubic, etc...)
pdf(file = 'HRF_signals_plots_multi_stim_holly_lmm.pdf', onefile = TRUE)
my_results_multi_lmm<-fTCD_glm_multi_lmm(path=paste0(getwd(),"/Holly_fTCD_data_run1"),order=order)
dev.off()
#stop()
#-----------------------------------------------------------------------------------------------------------------------#
#-----------------------------------------------------------------------------------------------------------------------#


mylong_resultsLMM<-gather(my_results_multi_lmm,key='param',value='beta',-c(ID,HRF,p_value_contrast))

names(mylong_resultsLMM)[2]<-"HRF"


ggplot(mylong_resultsLMM,aes(x=beta))+geom_density(fill='blue',alpha=0.5)+facet_wrap(~param,scales='free')+theme_bw()


