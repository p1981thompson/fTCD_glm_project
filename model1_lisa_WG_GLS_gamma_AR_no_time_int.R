#---------------------------------------------------------------------------------------------------------------#
# GLM model (without contrasts) for LISA data (Word Generation) - simple stimulus (stim1 and stim2) - both runs
#
# Model adapted from standard GLM to use GLS. Potential problem involving the induced autocorrelation from the heartbeat correction. Downsample the model to only use unique values (first value) from each block provides a better fit overall. Also, reorder the data so that the left and right signals are not sequential (i.e. left then right. We change the time variables so that left and right have the same range not 0-1, and 1-2, only 0-1).

# Original [NOT downsampled to 5Hz]
# Fixed single Gamma HRF


#--------------------------------------------------------------------------------------------------#

# Created by Paul Thompson and Zoe Woodhead - 17th Oct 2019
# Edited by Paul Thompson - 26th Nov 2019
# Edited by Paul Thompson - 30th June 2020
#--------------------------------------------------------------------------------------------------#

# install required packages for fitting the model.

#list_of_packages<-c("remotes","tidyverse","papaja","officer","fmri","knitr","utils","boot","ggpubr","psych","Rcpp","cladoRcpp","nlme","plm")
#new.packages <- list_of_packages[!(list_of_packages %in% installed.packages()[,"Package"])]
#if(length(new.packages))install.packages(new.packages,dependencies = TRUE)

#remotes::install_github("centerforopenscience/osfr")
#--------------------------------------------------------------------------------------------------#

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
library(plm)

#--------------------------------------------------------------------------------------------------#

#Load the data from Open Science Framework. 
#This can be done manually (got to: https://osf.io/5kq42/), or via the script below if user is comfortable with API tokens.

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
#---------------------------------------------------------------------------------------------------#

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
# PART 3:
#
#   - plot the correlagram
#

fTCD_glm_LISA_GLS<-function(path,order)
{
  # get all files names to be loaded in and preprocessed
  filename1<-list.files(path,pattern = '.exp')
  
  ## Set parameters
  samplingrate <- 25 # Sampling rate after downsampling. Raw data is 100Hz, we take 1 in every 4 samples
  heartratemax <- 125
  
  # set up data.frame to hold the outputted parameter estimates from the GLMs.
  glm.data<-data.frame(matrix(NA,nrow=length(filename1),ncol=(((order+5))+2)))
  names(glm.data)<-c('ID',paste0('param',(1:(order+5))),'HRF')
  
  #------------------------------------------------------------------------------------------------#
  #########################################
  # PART 1                                #
  #                                       #
  # Created by z.woodhead 30th July 2019  #
  # Edited  by z. woodhead 3rd Oct 2019   #
  #########################################
  #------------------------------------------------------------------------------------------------#

  for(j in 1:length(filename1))
  {
    print(filename1[j])
    
    ## Read in raw data
    
    myfile <- filename1[j]
    mydata<-read.table(paste0(path,"/",myfile), skip = 6,  header =FALSE, sep ='\t')
    
    wantcols = c(2,3,4,7) #sec, L, R,marker #select columns of interest to put in shortdat
    shortdat = data.frame(mydata[,wantcols])
    rawdata = filter(shortdat, row_number() %% 4 == 0) # downsample to 5 Hz by taking every 20th point
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
    
    #-----------------------------------------------------------------------------------------------#
    # Save processed file in csv format
    #mynewfile <- paste0(getwd(),"/Lisa_data/Chpt4_fTCD_WordGen_rawdata/",strsplit(myfile, '*.exp'), '_processed.csv')
    #write.csv(rawdata, mynewfile, row.names=F)
    
    #-----------------------------------------------------------------------------------------------#

    #########################################
    # PART 2                                #
    #                                       #
    # Created by P.Thompson 17th Oct 2019   #
    # Edited by P.Thompson 18th Oct 2019    #
    # Edited by P.Thompson 30th June 2020   #
    #########################################
    #-----------------------------------------------------------------------------------------------#
    rawdata2<-rawdata
    #-----------------------------------------------------------------------------------------------#
    
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
      
      .gammaHRF <- function(t, par = NULL) {
        th <- 0.242 * par[1]
        1/(th * factorial(3)) * (t/th)^3 * exp(-t/th)
      }
      
      
      par <- floor((durations[1]/28)*4)
      
      y <- .gammaHRF(0:(durations[1] * scale)/scale, par) 
      
      stimulus <-  convolve(stimulus, rev(y), type = "open")
      stimulus <- stimulus[unique((scale:scans)%/%(scale^2 * TR)) * scale^2 * TR]/(scale^2 * TR)
      stimulus <- stimulus - mean(stimulus)
      return(stimulus)
    }  
    
    #-----------------------------------------------------------------------------------------------#
    #-----------------------------------------------------------------------------------------------#
    
    # Create convolved stimulus function with HRF (applying the new fmri.stimulus.PT2 function above)
    
    gamma1 = fmri.stimulus.PT2(scans = dim(rawdata)[1], onsets = c(1,1+which(diff(rawdata$stim1_on)!=0)), durations = stim1_length_samples, TR = 1,scale=1)
    
    gamma2 = fmri.stimulus.PT2(scans = dim(rawdata)[1], onsets = c(1,1+which(diff(rawdata$stim2_on)!=0)), durations = stim2_length_samples, TR = 1,scale=1)
    
    #-----------------------------------------------------------------------------------------------# 
    # Binds all the stimuli into one matrix to be read into the fmri.design function. This converts the data into a design matrix and adds in the drift terms according to the order argument specified by the user.
    gamma = as.matrix(cbind(gamma1,gamma2))
    gamma = rbind(gamma,gamma)
    
    #-----------------------------------------------------------------------------------------------#
    
    # We create the design matrix and bind them together to give the same design matrix for each side (left and right), so that the main effect of side can be modelled appropriately.
    my_des<-fmri.design(gamma, order = order)
    
    # Add a dummy variable for side (signal). This is either 0 or 1 for left and right respectively.
    my_des<-cbind(my_des,rep(1:0,each=length(gamma1)))
    
    # Add a dummy variable for side (signal). This is either 0 or 1 for left and right respectively.
    my_des<-cbind(my_des,rep(1:0,each=length(gamma1)))
    
    # Add interaction variable for side (signal*stim1).
    my_des[,8]<-my_des[,8]*my_des[,1]
    
    #-----------------------------------------------------------------------------------------------#
    
    mydata<-data.frame(y=c(rawdata$heartbeatcorrected_L,rawdata$heartbeatcorrected_R),stim1=my_des[,1],stim2=my_des[,2],t=my_des[,4],signal=as.factor(my_des[,7]),stim1_signal=my_des[,8])
    
    mydata[(length(rawdata$heartbeatcorrected_L)+1):(length(rawdata$heartbeatcorrected_L)*2),4]<-mydata[(length(rawdata$heartbeatcorrected_L)+1):(length(rawdata$heartbeatcorrected_L)*2),4]-1
    
    #-----------------------------------------------------------------------------------------------#
    
    # Uses autocorrelated errors via gls model. (uses the 'nlme' package to achieve this error structure).
    
    # filter out replicates in the dependent variable relating the the heartbeat correction (articially induces autocorrelation if left in)
    mydata<-mydata %>% group_by(y) %>% filter(t==min(t))
    
    # set optimisation parameters 
    glsControl(optimMethod = "L-BFGS-B",maxIter = 100)
    
    # fit model using generalized least squares
    myfit <- nlme::gls(y~stim1+stim2+t+I(t^2)+I(t^3)+signal+stim1_signal,data=mydata,
                       correlation=corAR1(form=~1|signal))
    
    # Ensure class and coefficients are correctly labelled.
    names(myfit$coefficients)<-c("intercept","stim1","stim2","t","t_sqr","t_cub","signal","interaction")
    
    #-----------------------------------------------------------------------------------------------#
    
    # Extract the parameter estimates and record them for later use. Data stored in data.frame called 'glm.data'.
    glm.data[j,1] <- strsplit(basename(myfile),'[.]')[[1]][1]
    
    glm.data[j,(((order+5)*1)+2)] <- "gamma"
    
    glm.data[j,2:(((order+5)*1)+1)] <- myfit$coefficients
    
    #-----------------------------------------------------------------------------------------------#
    #setup plot data (wrangling data to work with plot)
    myplotdat<-data.frame(y=mydata$y,x=mydata$t,fitted=predict(myfit),Signal=mydata$signal)
    
    #plot the time series for each individual per signal.
    g3<-ggplot(myplotdat,aes(y=y,x=x,colour=Signal))+geom_point(aes(colour=Signal),alpha=0.4)+geom_line(aes(y=fitted))+theme_bw()
    
    # as we are fitting in a loop and printing to file, we need to use 'print' function with ggplot.
    print(g3)
    
    #Also plot the autocorrelation plot of the residuals
    acf(residuals(myfit,type="normalized"))
    
    #output data
    glm_data<-glm.data
  }
  
  return(glm_data)    
}

################################# END OF FUNCTION ###################################################


#--------------------------------------------------------------------------------------------------#
# RUN FUNCTION FOR ALL PARTICIPANT DATA FILES
#--------------------------------------------------------------------------------------------------#

#Set the order

order=3 #polynomial drift terms (2=quadratic, 3=cubic, etc...)
pdf(file = '/Users/paulthompson/Documents/fTCD_glm_results/Lisa_data/test models/HRF_signals_plots_LISA_WG_GLS_gamma_AR1.pdf', onefile = TRUE) #print plots to file.
my_results_LISA_WG_GLS_gamma_AR1<-fTCD_glm_LISA_GLS(path='/Users/paulthompson/Documents/fTCD_glm_results/Lisa_data/Chpt4_fTCD_WordGen_rawdata',order=order)
dev.off()

#---------------------------------------------------------------------------------------------------#
write.csv(my_results_LISA_WG_GLS_gamma_AR1,'/Users/paulthompson/Documents/fTCD_glm_results/Lisa_data/test models/my_results_LISA_WG_GLS_gamma_AR1.csv',row.names = FALSE)

my_results_LISA_WG_GLS_gamma_AR1<-my_results_LISA_WG_GLS_gamma_AR1[complete.cases(my_results_LISA_WG_GLS_gamma_AR1), ]
#---------------------------------------------------------------------------------------------------#

#Exclusions

exclude_id<-c(paste0(c('013','031','102','108','120','121','125','129','134','139','141','142'),'DAC1'))
# 
my_results_LISA_WG_GLS_gamma_AR1_ex <- my_results_LISA_WG_GLS_gamma_AR1[!my_results_LISA_WG_GLS_gamma_AR1$ID %in% exclude_id,]

#---------------------------------------------------------------------------------------------------#
# --------------------------------------------------------------------------------------------------#

#load LI based on old doppler analysis method

old_res<-read.csv("/Users/paulthompson/Documents/fTCD_glm_results/Lisa_data/WordGen_results.csv")

old_res<-old_res %>% rename(ID=Filename)

compare_results<-merge(my_results_LISA_WG_GLS_gamma_AR1_ex,old_res,by='ID',all.x = T)

fmri_data <- read.csv('/Users/paulthompson/Documents/fTCD_glm_results/Lisa_data/Chapter5_fMRI_data.csv')

# Identify factors
factor_variables <- c('group_cat', 'group_lat', 'sex', 'hand_self_report', 'hand_QHP_cat', 'hand_EHI_cat', 'Old_fTCD_wg_cat', 'Old_fTCD_pptt_cat', 'fTCD_wg_cat', 'fTCD_pptt_cat')

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

pdf(file = '/Users/paulthompson/Documents/fTCD_glm_results/Lisa_data/test models/HRF_signals_plots_LISA_WG_GLS_gamma_AR1_correlations.pdf')

psych::pairs.panels(compare_results2[,c('fMRI_diff_wg_frontal','fMRI_diff_wg_temporal','fMRI_diff_wg_MCA','LI','param8')])

dev.off()
