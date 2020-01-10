#---------------------------------------------------------------------------------------------------------------#
# GLM model (without contrasts) for A2 data (Sentence Generation) - simple stimulus (stim1 and stim2) - both runs
#
# Model adapted from standard GLM to have autocorrelated residuals according to Worsley et al. (2002) A general statistical analysis for fMRI data. 
# Downsampled to 5Hz.
# Fixed canonical HRF
# nlme::gls( AR(1)) model 

#---------------------------------------------------------------------------------------------------------------#

# Created by Paul Thompson and Zoe Woodhead - 17th Oct 2019
# Edited by Paul Thompson - 25th Nov 2019

# install required packages for fitting the model.

#list_of_packages<-c("remotes","tidyverse","papaja","officer","fmri","knitr","utils","boot","ggpubr","psych","Rcpp","cladoRcpp","nlme")
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
library(nlme)

#---------------------------------------------------------------------------------------------------------------#

#Load the data from Open Science Framework. This can be done manually (got to: https://osf.io/5kq42/), or via the script below if user is comfortable with API tokens.

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
    osf_download() %>% unzip('A2_SG_Data.zip')
}
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
# PART 2:
#
#   - runs the glm
#   - saves the parameter estimates to data.frame.
#

fTCD_glm_A2_AC<-function(path,order)
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
  
  for(j in (c(1:(length(filename1)/2)*2)-1)) # skips the second repeated measure (D1 only), only fits for first 35 ids.
  {
    print(filename1[j])
    
    ## Read in raw data
    
    myfile <- filename1[j]
    mydata<-read.table(paste0(path,"/",myfile), skip = 6,  header =FALSE, sep ='\t')
    
    wantcols = c(2,3,4,9) #sec, L, R,marker #select columns of interest to put in shortdat
    shortdat = data.frame(mydata[,wantcols])
    rawdata = filter(shortdat, row_number() %% 20 == 0) # downsample to 5 Hz by taking every 20th point
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
    stim1_delay_sec <- 3
    stim1_delay_samples <- stim1_delay_sec * samplingrate
    stim1_length_sec <- 14
    stim1_length_samples <- stim1_length_sec * samplingrate
    
    stim2_delay_sec <- 17
    stim2_delay_samples <- stim2_delay_sec * samplingrate
    stim2_length_sec <- 6
    stim2_length_samples <- stim2_length_sec * samplingrate
    
    rest_length_sec <- 10
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
    mynewfile <- paste0('/Volumes/PSYHOME/PSYRES/pthompson/DVMB/fTCD_glm_project',"/A2_SG_data/",strsplit(myfile, '*.exp'), '_processed.csv')
    #write.csv(rawdata, mynewfile, row.names=F)
    
    #---------------------------------------------------------------------------------------------------------------#
    
    #---------------------------------------------------------------------------------------------------------------#
    #########################################
    # PART 2                                #
    #                                       #
    # Created by P.Thompson 17th Oct 2019   #
    # Edited by P.Thompson 18th Oct 2019    #
    #########################################
    #---------------------------------------------------------------------------------------------------------------#
    
    #myseq<-seq(1,length(rawdata[,1]),by=25)
    rawdata2<-rawdata#[myseq,]
    
    
    #---------------------------------------------------------------------------------------------------------------#
    
    # Adapted 'fmri.stimulus' function from the R package 'fmri'. This is a condensed version that only gives option of the gamma HRF and convolves the HRF to the stimli specificied earlier in this script
    
    fmri.stimulus.PT2<- function(scans = dim(rawdata)[1], onsets = c(1,1+which(diff(rawdata$stim1_on)!=0)), durations = 375, TR = 1/25,scale=1)
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
      
      .canonicalHRF <- function(t, par = NULL) {
        ttpr <- par[1] * par[3]
        ttpu <- par[2] * par[4]
        (t/ttpr)^par[1] * exp(-(t - ttpr)/par[3]) - par[5] * 
          (t/ttpu)^par[2] * exp(-(t - ttpu)/par[4])
      }
      
      
      par <- c(6, 12, 0.9, 0.9, 0.35)
      
      y <- .canonicalHRF(0:(20 * scale)/scale, par)/2.885802
      
      stimulus <- convolve(stimulus, rev(y), type = "open")
      stimulus <- stimulus[unique((scale:scans)%/%(scale^2 * TR)) * scale^2 * TR]/(scale^2 * TR)
      stimulus <- stimulus - mean(stimulus)
      return(stimulus)
    }  
    
    
    #---------------------------------------------------------------------------------------------------------------#
    #---------------------------------------------------------------------------------------------------------------#
    # Create convolved stimulus function with HRF (applying the new fmri.stimulus.PT2 function above)
    
    gamma1 = fmri.stimulus.PT2(scans = dim(rawdata)[1], onsets = c(1,1+which(diff(rawdata$stim1_on)!=0)), durations = stim1_length_samples, TR = 1,scale=1)
    
    gamma2 = fmri.stimulus.PT2(scans = dim(rawdata)[1], onsets = c(1,1+which(diff(rawdata$stim2_on)!=0)), durations = stim2_length_samples, TR = 1,scale=1)
    
    #---------------------------------------------------------------------------------------------------------------# 
    # Binds all the stimuli into one matrix to be read into the fmri.design function. THis converts the data into a design matrix and adds in the drift terms according to the order argument specified by the user.
    gamma = as.matrix(cbind(gamma1,gamma2))
    
    # Binds all the stimuli into one matrix to be read into the fmri.design function. THis converts the data into a design matrix and adds in the drift terms according to the order argument specified by the user.
    gamma = rbind(gamma,gamma)
    
    #---------------------------------------------------------------------------------------------------------------#
    
    # We create the design matrix and bind them together to give the same design matrix for each side (left and right), so that the main effect of side can be modelled appropriately.
    my_des<-fmri.design(gamma, order = order)
    
    # Add a dummy variable for side (signal). This is either 0 or 1 for left and right respectively.
    my_des<-cbind(my_des,rep(1:0,each=length(gamma1)))
    
    # Add a dummy variable for side (signal). This is either 0 or 1 for left and right respectively.
    my_des<-cbind(my_des,rep(1:0,each=length(gamma1)))
    
    # Add interaction variable for side (signal*stim1).
    my_des[,8]<-my_des[,8]*my_des[,1]
    
    #---------------------------------------------------------------------------------------------------------------#
    
    mydata<-data.frame(y=c(rawdata$heartbeatcorrected_L,rawdata$heartbeatcorrected_R),stim1=my_des[,1],stim2=my_des[,2],t=my_des[,4],signal=as.factor(my_des[,7]),stim1_signal=my_des[,8])
    
    #---------------------------------------------------------------------------------------------------------------#
    # Uses autocorrelated errors via gls model. (uses the 'nlme' package to achieve this error structure).
    
    myfit <- gls(y~stim1+stim2+t+I(t^2)+I(t^3)+signal+stim1_signal,data=mydata,
                 correlation=corAR1(form=~t))
    
    print(names(myfit$coefficients))
    
    # Ensure class and coefficients are correctly labelled.
    names(myfit$coefficients)<-c("intercept","stim1","stim2","t","t_sqr","t_cub","signal","interaction")
    
    #---------------------------------------------------------------------------------------------------------------#
    
    # Extract the parameter estimates and record them for later use. Data stored in data.frame called 'glm.data'.
    glm.data[j,1] <- strsplit(basename(myfile),'[.]')[[1]][1]
    
    glm.data[j,(((order+5)*1)+2)] <- "gamma"
    
    glm.data[j,2:(((order+5)*1)+1)] <- myfit$coefficients
    
    #---------------------------------------------------------------------------------------------------------------#
    #setup data for plotting in ggplot
    pframe<-with(rawdata,expand.grid(t=seq(min(sec),max(sec),length=length(rawdata$heartbeatcorrected_L)),signal=c(0,1)))
    
    pframe<-data.frame(stim1=c(gamma1,gamma1),stim2=c(gamma2,gamma2),t=pframe[,1],t_sqr=(pframe[,1])^2,t_cub=(pframe[,1])^3,signal=pframe[,2],interaction=c(gamma1,gamma1)*pframe[,2])
    
    myplotdat<-data.frame(y=c(rawdata$heartbeatcorrected_L,rawdata$heartbeatcorrected_R),
                          x=c(rawdata$sec,rawdata$sec),
                          fitted=predict(myfit),Signal=rep(c("Left","Right"),each=length(rawdata$sec)))
    
    g3<-ggplot(myplotdat,aes(y=y,x=x,colour=Signal))+geom_point(colour='grey',alpha=0.5)+geom_line(aes(y=fitted))+theme_bw()
    
    # as we are fitting in a loop and printing to file, we need to use 'print' function with ggplot.
    print(g3)
    
    #output data
    glm_data<-glm.data
  }
  
  return(glm_data)    
}

################################# END OF FUNCTION #######################################################################


#-----------------------------------------------------------------------------------------------------------------------#
# RUN FUNCTION FOR ALL PARTICIPANT DATA FILES
#-----------------------------------------------------------------------------------------------------------------------#
#Set the order
order=3 #polynomial drift terms (2=quadratic, 3=cubic, etc...)
pdf(file = '/Volumes/PSYHOME/PSYRES/pthompson/DVMB/fTCD_glm_project/HRF_signals_plots_A2project_SG_AutoCor_Err_Downsampled5Hz_canonical.pdf', onefile = TRUE) #print plots to file.
my_results_A2_SG_AutoCor_Err_downsampled5Hz_canonical<-fTCD_glm_A2_AC(path=paste0('/Volumes/PSYHOME/PSYRES/pthompson/DVMB/fTCD_glm_project','/A2_SG_data'),order=order)
dev.off()

my_results_A2_SG_AutoCor_Err_downsampled5Hz_canonical<-my_results_A2_SG_AutoCor_Err_downsampled5Hz_canonical[complete.cases(my_results_A2_SG_AutoCor_Err_downsampled5Hz_canonical), ]

write.csv(my_results_A2_SG_AutoCor_Err_downsampled5Hz_canonical,'/Volumes/PSYHOME/PSYRES/pthompson/DVMB/fTCD_glm_project/downsample_at_5Hz/my_results_A2_SG_AutoCor_Err_downsampled5Hz_canonical.csv',row.names=FALSE)
#-----------------------------------------------------------------------------------------------------------------------#

exclude_id<-c(paste0('A2_',c('013','031'),'_D1'))
# 
my_results_A2_SG_AutoCor_Err_downsampled5Hz_canonical <- my_results_A2_SG_AutoCor_Err_downsampled5Hz_canonical[!my_results_A2_SG_AutoCor_Err_downsampled5Hz_canonical$ID %in% exclude_id,]
# 
# 
# 
# #-----------------------------------------------------------------------------------------------------------------------#
# --------------------------------------------------------------------------------------------------------------#
# 
my_results_A2_SG_AutoCor_Err_downsampled5Hz_canonical$ID <- substring(my_results_A2_SG_AutoCor_Err_downsampled5Hz_canonical$ID,4,6)
# 
# 
# #-----------------------------------------------------------------------------------------------------------------------#
# #load LI based on old doppler analysis method
# 
old_res_A2<-read.csv("/Volumes/PSYHOME/PSYRES/pthompson/DVMB/fTCD_glm_project/A2_SG_LI.csv")
# 
old_res_A2$ID<-sprintf('%0.3d', old_res_A2$ID)
# 
old_res_A2<-old_res_A2[,1:3]
# 
names(old_res_A2)<-c('ID','LI_SG1','LI_SG2')
# 
compare_results_A2_SG_AutoCor_Err_downsampled5Hz_canonical<-merge(my_results_A2_SG_AutoCor_Err_downsampled5Hz_canonical,old_res_A2,by='ID',all.x = T)
# 
# #-----------------------------------------------------------------------------------------------------------------------#
# 
# #-----------------------------------------------------------------------------------------------------------------------#
# 
# #Print correlation matrix plots to check association between the old LI and new glm-derived LI measures.
psych::pairs.panels(compare_results_A2_SG_AutoCor_Err_downsampled5Hz_canonical[,c('param8','LI_SG1')],cex.cor=1)
# 
# #-----------------------------------------------------------------------------------------------------------------------#
