# Comparing centWave to milliWave

library(xcms)

source("xcms/milliWave_source.R")

load("xcms/raw_data")

peakwidth <- c(20, 80) # Minimum and maximum width of a few random peaks, in seconds
mz_span <- c(0.0005) # Maximum spread of m/z values across a well-defined peak, plus some buffer
ppm <- ceiling((mz_span*1000000)/132)
cwp <- CentWaveParam(peakwidth = peakwidth, ppm = ppm, 
                     snthresh = 1, prefilter = c(3,10000))

# xdata <- findChromPeaks(raw_data, param = cwp) #Takes about 25 mins in serial on laptop
# save(xdata, file = "xcms/xdata")
load("xcms/xdata")


mdata <- will_findChromPeaks_milliWave(raw_data, param = cwp)
save(mdata, file = "xcms/mdata/milliWave_code")
load("xcms/mdata/milliWave_code")
