library(fmri)
library(tidyverse)
library(splines)
library(boot)

#--------------------------------------------------------------------------------------#

# We use the framework from the fmri based analysis, but substitute the Lindquist et al 2009 HRF in to the function. The data passed through is actually fTCD data, so parameters for 'scans', 'onsets', 'durations', to something used in doppler which were the sampling rate, time between sampling and duration of each epoch.

#' @article{LINDQUIST2009S187,
#'   title = "Modeling the hemodynamic response function in fMRI: Efficiency, bias and mis-modeling",
#'   journal = "NeuroImage",
#'   volume = "45",
#'   number = "1, Supplement 1",
#'   pages = "S187 - S198",
#'   year = "2009",
#'   note = "Mathematics in Brain Imaging",
#'   issn = "1053-8119",
#'   doi = "https://doi.org/10.1016/j.neuroimage.2008.10.065",
#'   url = "http://www.sciencedirect.com/science/article/pii/S1053811908012056",
#'   author = "Martin A. Lindquist and Ji Meng Loh and Lauren Y. Atlas and Tor D. Wager",
#' }


# https://www.sciencedirect.com/science/article/pii/S1053811908012056?via%3Dihub 


# citation(package = "fmri")
# Karsten Tabelow, Joerg Polzehl (2011). Statistical Parametric Maps for Functional MRI Experiments in R: The Package fmri. Journal of Statistical
# Software, 44(11), 1-21. URL http://www.jstatsoft.org/v44/i11/.


Lindquist_HRF <- function(t, par = NULL) {
  require(splines)
  fm1 <- lm(par ~ bs(t, df = 3))
  
  ## example of safe prediction
  ht <- seq(0, max(t), length.out = max(t))
  return(predict(fm1, data.frame(t = ht)))
  
}

# -------------------------------------------------------------------------------------#



# Functions from the fmri package, with parameters adapted for fTCD input.

#Boxcar

fTCD.bold.resp1 <- fmri.stimulus(scans = dim(data1)[1], onsets = c(1,1+which(diff(data1$stim_on)!=0))[seq(2, length(c(1,1+which(diff(data1$stim_on)!=0))), by = 2)], durations = 500, TR = 1/25,type="boxcar",scale=1)

fTCD.bold.resp2 <- fmri.stimulus(scans = dim(data1)[1], onsets = c(1,1+which(diff(data1$stim_on)!=0))[seq(2, length(c(1,1+which(diff(data1$stim_on)!=0))), by = 2)], durations = 500, TR = 1/25,type="canonical",scale=1)


my_des<-fmri.design(fTCD.bold.resp1)

par(mfrow=c(2, 2))
for (i in 1:4) plot(my_des[, i], type="l")
par(mfrow=c(1, 1))

myfit<-glm.fit(x=my_des,y=data1$heartbeatcorrected_L[c(seq(from=1, to=length(data1[,1]), by=25)-1)[-1]],family=gaussian())

myfit$coefficients

class(myfit) <- c(myfit$class, c("glm", "lm"))

myfit.diag<-glm.diag(myfit)

glm.diag.plots(myfit, myfit.diag)
# -------------------------------------------------------------------------------------#

#canonical

fTCD.bold.resp2 <- fmri.stimulus(scans = dim(data1)[1], onsets = c(1,1+which(diff(data1$stim_on)!=0))[seq(2, length(c(1,1+which(diff(data1$stim_on)!=0))), by = 2)], durations = 500, TR = 1/25,type="canonical",scale=1)


my_des2<-fmri.design(fTCD.bold.resp2)


myfit2<-glm.fit(x=my_des2,y=data1$heartbeatcorrected_L[c(seq(from=1, to=length(data1[,1]), by=25)-1)[-1]])

myfit2$coefficients



class(myfit2) <- c(myfit2$class, c("glm", "lm"))

myfit.diag2<-glm.diag(myfit2)

glm.diag.plots(myfit2, myfit.diag2)

# -------------------------------------------------------------------------------------#

#gamma

fTCD.bold.resp3 <- fmri.stimulus(scans = dim(data1)[1], onsets = c(1,1+which(diff(data1$stim_on)!=0))[seq(2, length(c(1,1+which(diff(data1$stim_on)!=0))), by = 2)], durations = 500, TR = 1/25,type="gamma",scale=1)


my_des3<-fmri.design(fTCD.bold.resp3)


myfit3<-glm.fit(x=my_des3,y=data1$heartbeatcorrected_L[c(seq(from=1, to=length(data1[,1]), by=25)-1)[-1]])

myfit3$coefficients

class(myfit3) <- c(myfit3$class, c("glm", "lm"))

myfit.diag3<-glm.diag(myfit3)

glm.diag.plots(myfit3, myfit.diag3)

# -------------------------------------------------------------------------------------#
modelComp<-data.frame(HRF=c("Boxcar","Canonical","Gamma"),AIC=c(myfit$aic,myfit2$aic,myfit3$aic))
modelComp

# -------------------------------------------------------------------------------------#
# plot glm on real data.

myplotdat<-data.frame(y=data1$heartbeatcorrected_L[c(seq(from=1, to=length(data1[,1]), by=25)-1)[-1]],
                      x=data1$sec[c(seq(from=1, to=length(data1[,1]), by=25)-1)[-1]],
                      fitted=myfit3$fitted.values)

ggplot(myplotdat,aes(y=y,x=x))+geom_point()+geom_line(aes(y=fitted),color="red")+theme_bw()



# -------------------------------------------------------------------------------------#


