library(fmri)
library(tidyverse)
library(splines)

Lindquist_HRF <- function(t, par = NULL) {
  require(splines)
  fm1 <- lm(par ~ bs(t, df = 3))
  
  ## example of safe prediction
  ht <- seq(0, max(t), length.out = max(t))
  return(predict(fm1, data.frame(t = ht)))
  
}


fTCD.stimulus.PT <- function (scans = 1, onsets = c(1), durations = c(1), TR = 2, 
                              times = FALSE, sliceorder = NULL, type = c("canonical", "gamma", 
                                                                         "boxcar", "user"), par = NULL, scale = 1, hrf = NULL, 
                              verbose = FALSE) 
{
  type <- match.arg(type)
  if ((type == "user") && (class(hrf) != "function")) 
    stop("HRF type is user, but specified hrf is not a function!")
  if (!times) {
    onsets <- onsets * TR
    durations <- durations * TR
  }
  slicetiming <- !is.null(sliceorder)
  if (slicetiming) {
    nslices <- length(sliceorder)
    scale <- max(scale, nslices)
    slicetimes <- (1:nslices - 1)[sliceorder]/TR * scale
  }
  onsets <- onsets * scale
  durations <- durations * scale
  scans <- scans * TR * scale
  TR <- TR/scale
  slicetiming <- !is.null(sliceorder)
  if (slicetiming) {
    nslices <- length(sliceorder)
    slicetimes <- ceiling((1:nslices - 1)[sliceorder]/nslices * 
                            scale)
  }
  if (type == "user") 
    shrf <- sum(hrf(0:(ceiling(scans) - 1)/scale))
  no <- length(onsets)
  if (length(durations) == 1) {
    durations <- rep(durations, no)
  }
  else if (length(durations) != no) {
    stop("Length of duration vector does not match the number of onsets!")
  }
  if (slicetiming) {
    stimulus <- matrix(0, ceiling(scans), nslices)
    for (j in 1:nslices) for (i in 1:no) stimulus[pmax(1, 
                                                       onsets[i]:(onsets[i] + durations[i] - 1) - slicetimes[j]), 
                                                  j] <- 1
  }
  else {
    stimulus <- rep(0, ceiling(scans))
    for (i in 1:no) stimulus[onsets[i]:(onsets[i] + durations[i] - 
                                          1)] <- 1
  }
  .canonicalHRF <- function(t, par = NULL) {
    ttpr <- par[1] * par[3]
    ttpu <- par[2] * par[4]
    (t/ttpr)^par[1] * exp(-(t - ttpr)/par[3]) - par[5] * 
      (t/ttpu)^par[2] * exp(-(t - ttpu)/par[4])
  }
  .gammaHRF <- function(t, par = NULL) {
    th <- 0.242 * par[1]
    1/(th * factorial(3)) * (t/th)^3 * exp(-t/th)
  }
  if (type == "canonical") {
    if (is.null(par)) 
      par <- c(6, 12, 0.9, 0.9, 0.35)
    if (!is.numeric(par[1:5]) || any(is.na(par[1:5]))) {
      warning("parameter vector c(", paste(par, collapse = ", "), 
              ") for canonical HRF is not numeric or has unsufficient length (<5)!\nUsing default parameters!", 
              paste(par <- c(6, 12, 0.9, 0.9, 0.35), collapse = ", "))
    }
  }
  else if (type == "gamma") {
    if (is.null(par)) 
      par <- 4
    if (!is.numeric(par[1])) {
      warning("parameter ", par[1], " for gamma HRF is not numeric!\nUsing default parameter!", 
              par <- 4)
    }
  }
  if (verbose) 
    cat("fmriStimulus: Using", type, "HRF for stimulus creation\n")
  y <- switch(type, canonical = .canonicalHRF(0:(20 * scale)/scale, 
                                              par)/2.885802, gamma = .gammaHRF(0:(28 * scale)/scale, 
                                                                               par), boxcar = scale, user = hrf(0:(ceiling(scans) - 
                                                                                                                     1)/scale)/shrf,par)
  if (slicetiming) {
    for (j in 1:nslices) {
      stimulus[, j] <- convolve(stimulus[, j], rev(y), 
                                type = "open")[1:dim(stimulus)[1]]
    }
    stimulus <- stimulus[unique((scale:scans)%/%(scale^2 * 
                                                   TR)) * scale^2 * TR, ]/(scale^2 * TR)
  }
  else {
    stimulus <- convolve(stimulus, rev(y), type = "open")
    stimulus <- stimulus[unique((scale:scans)%/%(scale^2 * 
                                                   TR)) * scale^2 * TR]/(scale^2 * TR)
  }
  if (slicetiming) 
    sweep(stimulus, 2, apply(stimulus, 2, mean), "-")
  else stimulus - mean(stimulus)
}


#--------------------------------------------------------------------------------------#



fTCD.stimulus.PT(scans = 1, onsets = c(1), durations = c(1), TR = 2,type="user",hrf=Lindquist_HRF,par='path/datafile location')


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
