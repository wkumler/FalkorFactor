# NoiseAlternateFunctions.R
# Contains alternatives to xcms:::getLocalNoiseEstimate

getLocalNoise_maxpeaks <- function(d, peakwidth){
  noise <- head(sort(d), -max(peakwidth))
  return(c(mean(noise), sd(noise)))
}

getLocalNoise_sd <- function(d){
  noise <- d[d<(mean(d)+sd(d))]
  return(c(mean(noise), sd(noise)))
}

getLocalNoise_IQR <- function(d){
  noise <- d[d<(median(d)+IQR(d))]
  return(c(mean(noise), sd(noise)))
}