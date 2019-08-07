


library(tidyverse)
library(forecast)

#============================================================================================================================#
#read in the data file

data1<-read_csv("013DAC1_processed.csv")
data1<-as.data.frame(data1)

data1$stim_on<-car::recode(data1$stim_on,"0=100;1=110")

data2<-read_csv("113DAC1_processed.csv")
data2<-as.data.frame(data2)

data2$stim_on<-car::recode(data2$stim_on,"0=100;1=110")

#library(signal)

library(eegkit)

data1$filtered_L = eegfilter(data1$heartbeatcorrected_L, Fs = 25, lower = 1/100, upper = 5, method = "butter", order = 4)
data2$filtered_L = eegfilter(data2$heartbeatcorrected_L, Fs = 25, lower = 1/100, upper = 5, method = "butter", order = 4)
# bf <- butter(2,1/1200*2, type="low")
# data1$filtered_L = signal:::filter(bf, data1$heartbeatcorrected_L)
# data2$filtered_L = signal:::filter(bf, data2$heartbeatcorrected_L)
#============================================================================================================================#

long_ts1<-gather(data1,key='signal',value = 'y',-sec)
long_ts2<-gather(data2,key='signal',value = 'y',-sec)

long_ts1$id<-rep("013DAC1",length(long_ts1[,1]))
long_ts2$id<-rep("113DAC1",length(long_ts2[,1]))

long_ts<-rbind(long_ts1,long_ts2)

long_ts2<-long_ts %>% dplyr::filter(.,signal=='heartbeatcorrected_L')
long_ts2a<-long_ts %>% dplyr::filter(.,signal=='stim_on')

ggplot(long_ts2, aes(x = sec, y = y)) + 
  geom_line(data=long_ts2a,aes(x=sec,y=y,linetype=id),color='black',alpha=0.7)+
  geom_line(aes(color=id))+
  theme_minimal()+theme(legend.position='top') 

#============================================================================================================================#

long_ts2<-long_ts %>% dplyr::filter(.,signal=='heartbeatcorrected_L'|signal=='stim_on')

fTCD.bold.resp2 <- fmri.stimulus(scans = dim(data)[1], onsets = c(1,1+which(diff(data1$stim_on)!=0))[seq(2, length(c(1,1+which(diff(data1$stim_on)!=0))), by = 2)], durations = 500, TR = 1/25,type="boxcar",scale=1)

ggplot(long_ts2, aes(x = sec, y = y)) + 
  geom_line(aes(color=signal))+facet_grid(id~.)+
  theme_minimal()+theme(legend.position='top') + geom_line(data=data.frame(y=fTCD.bold.resp2+100,x=seq(1:1239)),aes(y=y,x=x),color="black") 

#============================================================================================================================#

#filter time series signal to remove some noise.

long_ts3<-long_ts %>% dplyr::filter(.,signal=='filtered_L')
long_ts3a<-long_ts %>% dplyr::filter(.,signal=='stim_on')

ggplot(long_ts3, aes(x = sec, y = y)) + 
  geom_line(data=long_ts3a,aes(x=sec,y=y,linetype=id),color='black',alpha=0.7)+
  geom_line(aes(color=id))+
  theme_minimal()+theme(legend.position='top') 

#============================================================================================================================#
#============================================================================================================================#
#============================================================================================================================#
#============================================================================================================================#
#============================================================================================================================#
