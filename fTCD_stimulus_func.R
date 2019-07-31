library(fmri)
library(tidyverse)
library(splines)

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

# Functions from the fmri package, with parameters adapted for fTCD input.

fTCD.bold.resp <- fmri.stimulus(scans = 1, onsets = c(1), durations = c(1), TR = 2,type="user",hrf=Lindquist_HRF,par='path/datafile location',scale=1)



# fmri.stimulus(scans = 1, onsets = c(1), durations = c(1), TR = 2,
#               times = FALSE, sliceorder = NULL,
#               type = c("canonical", "gamma", "boxcar", "user"),
#               par = NULL, scale = 10, hrf = NULL, verbose = FALSE)

# -------------------------------------------------------------------------------------#
# scans	- number of scans
# 
# onsets - vector of onset times (in scans)
# 
# durations	- vector of duration of ON stimulus in scans or seconds (if !is.null(times))
# 
# TR - time between scans in seconds (TR)
# 
# times -	onset times in seconds. If present onsets arguments is ignored.
# 
# sliceorder - order of slice acquisition. If provided separate expected bold responses are calculated for the slices taking slice acquisition times into account. Default: no slice timing.
# 
# type - One of "canonical", "gamma", "boxcar", "user"
# 
# par	- Possible parameters to the HRF.
# 
# scale	- Temporal undersampling factor
# 
# hrf	- If type is "user" this should be a function evaluating the hemodynamic response function
# 
# verbose	- Report more if TRUE

# -------------------------------------------------------------------------------------#
