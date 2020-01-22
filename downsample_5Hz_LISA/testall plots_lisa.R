#---------------------------------------------------------------------------------------------------------------#
# GLM model (without contrasts) for A2 data (Sentence Generation) - simple stimulus (stim1 and stim2) - both runs
#---------------------------------------------------------------------------------------------------------------#
# We use generalised least squares to allow for the temporal dependence in the time series data. Each individual has
# two time series (left and right) measuring the response to stimulus in the middle cerebral arteries. The Hemodynamic 
# response has been calibrated using data from test subjects and we find corresponding evidence that the canonical 
# function (double gamma) gives a good fit for single response stimuli, but we can adapt this for sustained response 
# tasks and estimate the hemodynamic response function using a set of basis functions.
#
#---------------------------------------------------------------------------------------------------------------#

# Created by Paul Thompson and Zoe Woodhead - 17th Oct 2019

# install required packages for fitting the model.

#list_of_packages<-c("remotes","tidyverse","papaja","officer","fmri","knitr","utils","boot","ggpubr","psych","Rcpp","cladoRcpp","RcppArmadillo","doParallel","foreach")
#new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
#if(length(new.packages))install.packages(new.packages,dependencies = TRUE)

#remotes::install_github("centerforopenscience/osfr")


## Packages
library(osfr)
library(utils)
require(dplyr)
require(tidyverse)
require(boot)
require(fmri)
require(ggpubr)
library(psych)
library(cladoRcpp) # turbo charges the 'convolve' function used to convolve the stimulus and HRF (super slow previously). Based on the function 'fmri.stimulus' from the fmri package.
library(Rcpp)
library(RcppArmadillo)
library(doParallel)
library(foreach)

#---------------------------------------------------------------------------------------------------------------#

# This is the main function to run the analysis. The function does the following in order:
#
# PART 1:
#
#   Script takes a raw .exp datafile and preprocesses it ready for GLM analysis:
#   - It creates a box car function showing when the task was ON or OFF
#   - It normalises the fTCD signal to a mean of 100
#   - It performs heart beat integration
#   - It saves the processed data into a .csv file
#
# PART 1A:
#   Script epochs the processed raw data to fit the HRF function
#   - It epochs the continuous, processed data (rawdata) into trials (myepoched)
#   - It baseline corrects the epoched data
#   - It averages over all trials (myepoched_average)
#   - It plots the averaged data
#
# PART 2:
#
#   - runs the glm using the processed (continuous) data from part 1, and the HRF parameters from part 1A
#   - saves the parameter estimates to data.frame.
#

fTCD_glm_lisa_rcpp<-function(path,order)
{
  # get all files names to be loaded in and preprocessed
  filename1<-list.files(path,pattern = '.exp')
  
  ## Set parameters
  samplingrate <- 5 # Sampling rate after downsampling. Raw data is 100Hz, we take 1 in every 4 samples
  heartratemax <- 125
  
  # set up data.frame to hold the outputted parameter estimates from the GLMs.
  glm.data<-data.frame(matrix(NA,nrow=length(filename1),ncol=(((order+5))+2)))
  names(glm.data)<-c('ID',paste0('param',(1:(order+5))),'HRF')
  
  #---------------------------------------------------------------------------------------------------------------#
  #########################################
  # PART 1                                #
  #                                       #
  # Created by z.woodhead 30th July 2019  #
  # Edited  by z. woodhead 3rd Oct 2019   #
  #########################################
  #---------------------------------------------------------------------------------------------------------------#
  
  for(j in 1:length(filename1))
  {
    
    # cl <- parallel::makeForkCluster(2)
    # doParallel::registerDoParallel(cl)
    # foreach(j = 1:length(filename1), .combine = 'c') %dopar% {
    #   
    print(filename1[j])
    
    ## Read in raw data
    
    myfile <- filename1[j]
    prep_me<-1
    
    if(prep_me==1)
    {
      mydata<-read.table(paste0(path,"/",myfile), skip = 6,  header =FALSE, sep ='\t')
      
      wantcols = c(2,3,4,7) #sec, L, R,marker #select columns of interest to put in shortdat
      shortdat = data.frame(mydata[,wantcols])
      rawdata = filter(shortdat, row_number() %% 20 == 0) # downsample to 25 Hz by taking every 4th point
      allpts = nrow(rawdata) # total N points in long file
      rawdata[,1] = (seq(from=1,to=allpts*20,by=20)-1)/100 #create 1st column which is time in seconds from start
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
      
      #---------------------------------------------------------------------------------------------------------------#
      # Save processed file in csv format
      #mynewfile <- paste0(getwd(),"/A2_SG_data/",strsplit(myfile, '*.exp'), '_processed.csv')
      #write.csv(rawdata, mynewfile, row.names=F)
    }
    
    #rawdata <- read.csv(paste0("/Volumes/PSYHOME/PSYRES/pthompson/DVMB/fTCD_glm_project/A2_SG_data/",strsplit(myfile, '*.exp'), '_processed.csv'))
    
    #---------------------------------------------------------------------------------------------------------------#
    
    #---------------------------------------------------------------------------------------------------------------#
    #########################################
    # PART 1A                             #
    #                                     #
    # Created by Z.Woodhead 6th Nov 2019  #
    #########################################
    #---------------------------------------------------------------------------------------------------------------#
    
    # Epoch timings
    epochstart_time   <- -12
    epochend_time     <- 30
    epochstart_index  <- epochstart_time * samplingrate
    epochend_index    <- epochend_time * samplingrate
    basestart_time    <- -10 # baseline start
    baseend_time      <- 0 # baseline end
    basestart_index   <- basestart_time * samplingrate
    baseend_index    <- baseend_time * samplingrate
    
    # myepoched will be the full epoched trial
    myepoched <- array(0, dim=c(norigmarkers,epochend_index - epochstart_index + 1, 2)) # Set up an empty matrix
    
    for(mym in 1:norigmarkers) # for trials
    { 
      index1 = origmarkerlist[mym] + epochstart_index # index1 is index of the timepoint at the start of the epoch
      index2 = origmarkerlist[mym] + epochend_index # index2 is the index of the timepoint at the end of the epoch
      
      # If recording started late, the start of the epoch for trial 1 will be beyond the recorded range. 
      # If this doesn't affect the baseline period (ie, results will be unaffected), then replace with mean
      if (index1 < 0 & origmarkerlist[mym] + basestart_index > 0){
        cat("Recording started late. Padding start with zeros", "\n")
        replacement_mean_left = mean(rawdata[0 : index2, 2]) # Left hemisphere mean
        replacement_mean_right = mean(rawdata[0 : index2, 3]) # Right hemisphere mean
        # The epoched data is the heartbeat corrected data (columns 9 and 10 from rawdata)
        myepoched[mym, ,1] = c(rep(replacement_mean_left,index1*-1+1),rawdata[0:index2,9])
        myepoched[mym, ,2] = c(rep(replacement_mean_right,index1*-1+1),rawdata[0:index2,10])
      }
      
      if (index1 > 1){
        myepoched[mym,,1]=rawdata[index1:index2,9] #L side
        myepoched[mym,,2]=rawdata[index1:index2,10] #R side
      }
    }
    
    # Baseline correction
    basepoints=(basestart_index-epochstart_index):(baseend_index-epochstart_index) #all baseline points within epoch
    
    for (mym in 1:norigmarkers)
    {basemeanL=mean(myepoched[mym,basepoints,1]) #last dim is 3, which is HB corrected
    basemeanR=mean(myepoched[mym,basepoints,2])
    myepoched[mym,,1]=100+myepoched[mym,,1]-basemeanL #last dim 4 is HB and baseline
    myepoched[mym,,2]=100+myepoched[mym,,2]-basemeanR
    }
    
    # Average over trials
    ntime <- dim(myepoched)[2]
    myepoched_average <- data.frame(
      "Lmean" <- rep(1, ntime),
      "Rmean" <- rep(1, ntime),
      "LRdiff" <- rep(1, ntime))
    
    myepoched_average$Lmean <- apply(myepoched[ , , 1], c(2), mean)
    myepoched_average$Rmean <- apply(myepoched[ , , 2], c(2), mean)
    myepoched_average$LRdiff <- myepoched_average$Lmean - myepoched_average$Rmean
    myepoched_average$time<-seq(from=epochstart_time, to=epochend_time, by=.2)
    
    
    #---------------------------------------------------------------------------------------------------------------#
    #########################################
    # PART 2                                #
    #                                       #
    # Created by P.Thompson 17th Oct 2019  #
    #########################################
    #---------------------------------------------------------------------------------------------------------------#
    
    #myseq<-seq(1,length(rawdata[,1]),by=25)
    rawdata2<-rawdata#[myseq,]
    
    
    #---------------------------------------------------------------------------------------------------------------# 
    
    # Adapted 'fmri.stimulus' function from the R package 'fmri'. This is a condensed version that only gives option of the canonical HRF and estimates the individuals exact parameterisation of the canonical form. This then convolves the HRF to the stimli specificied earlier in this script
    
    
    fmri.stimulus.PT3 <- function(scans = dim(rawdata)[1], onsets = c(1,1+which(diff(rawdata$stim1_on)!=0)), durations = 375, TR = 1/25,scale=1,epoch_data=myepoched_average,startend=startend)
    {
      
      mylongdata<-gather(epoch_data[,c('time','Lmean','Rmean')],"signal","Mean_blood_flow",Lmean:Rmean)
      
      #---------------------------------------------------------------------------#
      # canonical function with parameter setup from:  https://doi.org/10.1016/j.neuroimage.2014.02.018
      
      mycanonicalHRF <- function(t, peak1,peak2,fwhm1,fwhm2,ratio,amp) {
        
        A <- 8*log(2)
        
        par1 <- A*((peak1^2)/(fwhm1^2))
        
        par2 <- (A*(peak2^2))/(fwhm2^2)
        
        par3 <- (fwhm1^2)/(A*peak1)
        
        par4 <- (fwhm2^2)/(A*peak2)
        
        ttpr <- par1 * par3
        ttpu <- par2 * par4
        
        resp <-100 + amp*((t/ttpr)^par1 * exp(-(t - ttpr)/par3) - ratio * (t/ttpu)^par2 * exp(-(t - ttpu)/par4))
        
        return(resp)
      }
      
      #---------------------------------------------------------------------------#
      
      # Use non-linear regression to fit the double gamma canonical function.
      
      library(minpack.lm)
      
      # create mean value between the Left and Right signals, so that the HRF is not over estimated.
      epoch_data$mean<-rowMeans(epoch_data[,c('Lmean','Rmean')])
      
      #manuallyset according to the stimuli start and end. Follwoing Zoe.W doppler script.
      timewindow1 <- startend[1]
      timewindow2 <- startend[2]
      
      #Automate the selection of starting values
      startamp <- max(epoch_data[epoch_data$time>=timewindow1 & epoch_data$time<=timewindow2,'mean'])-100
      startpeak1 <- epoch_data[epoch_data$time>=timewindow1 & epoch_data$time<=timewindow2,'time'][which(epoch_data[epoch_data$time>=timewindow1 & epoch_data$time<=timewindow2,'mean']==max(epoch_data[epoch_data$time>=timewindow1 & epoch_data$time<=timewindow2,'mean']))[1]]
      startpeak2 <- epoch_data[epoch_data$time>=timewindow1 & epoch_data$time<=timewindow2,'time'][which(epoch_data[epoch_data$time>=timewindow1 & epoch_data$time<=timewindow2,'mean']==min(epoch_data[epoch_data$time>=timewindow1 & epoch_data$time<=timewindow2,'mean']))[1]]
      
      if(startpeak2<=startpeak1+1){startpeak2<-startpeak1+2}
      
      startratio<-(min(epoch_data[epoch_data$time>=timewindow1 & epoch_data$time<=timewindow2,'mean'])[1])/(max(epoch_data[epoch_data$time>=timewindow1 & epoch_data$time<=timewindow2,'mean'])[1])
      
      
      
      # Fits a nonlinear modelling according to the model precedent set out in Proulx et al.(2014)
      #myHRFfit <- nlsLM(Mean_blood_flow~mycanonicalHRF(t=time,peak1,peak2,fwhm1,fwhm2,ratio,amp),data=mylongdata[mylongdata$time>=timewindow1 & mylongdata$time<=timewindow2,c(1,3)],start=list(peak1 = startpeak1, peak2 = startpeak2, fwhm1 = 5, fwhm2 = 5, ratio = startratio,amp=startamp),algorithm = "LM", lower=c(.01,.01,.01,.01,0.1,0.01),upper=c(50,50,10,10,2,10))
      
      print(ggplot(mylongdata,aes(y=Mean_blood_flow,x=time))+geom_line(aes(colour=signal))+geom_hline(yintercept=100,alpha=0.5)+geom_vline(xintercept = c(0,5,20),linetype='dashed',alpha=0.5)+theme_bw()+ggtitle(filename1[j]))#+ geom_line(data=data.frame(y=predict(myHRFfit,time=seq(timewindow1,timewindow2,by=0.1)),x=rep(seq(timewindow1,timewindow2,by=0.1),2)),aes(y=y,x=x),colour="purple"))
      
      #---------------------------------------------------------------------------#
      
    }  
    gamma1 = fmri.stimulus.PT3(scans = dim(rawdata)[1], onsets = c(1,1+which(diff(rawdata$stim1_on)!=0)), durations = stim1_length_samples, TR = 1,scale=1,epoch_data = myepoched_average,startend<-c(5,20))
    
    
  }
  
}

################################# END OF FUNCTION #######################################################################


#-----------------------------------------------------------------------------------------------------------------------#
# RUN FUNCTION FOR ALL PARTICIPANT DATA FILES
#-----------------------------------------------------------------------------------------------------------------------#
#Set the order
order=3 #polynomial drift terms (2=quadratic, 3=cubic, etc...)
pdf(file = 'HRF_signals_plots_LISA_WG_rcpp_downsamp5hz_TESTALL.pdf', onefile = TRUE) #print plots to file.
#start_time = Sys.time()
fTCD_glm_lisa_rcpp(path="/Volumes/PSYHOME/PSYRES/pthompson/DVMB/fTCD_glm_project/Chpt4_fTCD_WordGen_rawdata",order=order)
#end_time = Sys.time()

#end_time - start_time

dev.off()

#-----------------------------------------------------------------------------------------------------------------------#

#Exclusions

