#########################################################################################
# fTCD preprocessing for GLM analysis
# 
# Created by z.woodhead 30th July 2019 
# edited by P.Thompson 5th August 2019 - allow files to be read differently
#########################################################################################
# This script takes a raw .exp datafile and preprocesses it ready for GLM analysis:
#   - It creates a box car function showing when the task was ON or OFF
#   - It normalises the fTCD signal to a mean of 100
#   - It performs heart beat integration
#   - It saves the processed data into a .csv file

## Packages
require(dplyr)

## Set parameters
samplingrate <- 25 # Sampling rate after downsampling. Raw data is 100Hz, we take 1 in every 4 samples
heartratemax <- 125

## Read in raw data

myfile <- file.choose()
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
stim_delay_sec <- 5
stim_delay_samples <- stim_delay_sec * samplingrate
stim_length_sec <- 20
stim_length_samples <- stim_length_sec * samplingrate
rest_length_sec <- 30
rest_length_samples <- rest_length_sec * samplingrate

rawdata$stim_on <- 0
for (m in 1:norigmarkers){
  rawdata$stim_on[(origmarkerlist[m]+stim_delay_samples):(origmarkerlist[m]+stim_delay_samples+stim_length_samples)] <- 1
}

## Check markers look good - there should be 23 of them, and the time should run to over 1150 seconds (23 trials * 50 seconds)
plot(rawdata$sec, rawdata$stim_on, type='l')
cat("Press [enter] to continue")
line <- readline()


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

# Check data looks good
plot(rawdata$sec[1:400], rawdata$normal_L[1:400], type='l')
lines(rawdata$sec[1:400], rawdata$heartbeatcorrected_L[1:400], col='red')
cat("Press [enter] to continue")
line <- readline()

#----------------------------------------------------------
# Save processed file in csv format
mynewfile <- paste0(strsplit(myfile, '*.exp'), '_processed.csv')
write.csv(rawdata, mynewfile, row.names=F)
