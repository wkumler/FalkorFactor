# Comparing centWave to milliWave


# Setup things ----
setwd("C:/Users/willi/Documents/UW/FalkorFactor")
library(xcms)
source("xcms/milliWave_code/milliWave_source.R")
load("xcms/raw_data")

peakwidth <- c(20, 80) # Minimum and maximum width of a few random peaks, in seconds
mz_span <- c(0.0005) # Maximum spread of m/z values across a well-defined peak, plus some buffer
ppm <- ceiling((mz_span*1000000)/132)
cwp <- CentWaveParam(peakwidth = peakwidth, ppm = ppm, snthresh = 1, prefilter = c(3,10000))



# Compare code ----
sink("xcms/milliWave_code/centWaveCode.txt")
findMethods("findChromPeaks", classes = "OnDiskMSnExp#CentWaveParam")[[1]]
xcms:::findChromPeaks_OnDiskMSnExp
xcms:::findChromPeaks_Spectrum_list
do_findChromPeaks_centWave
xcms:::.centWave_new
sink()

sink("xcms/milliWave_code/milliWaveCode.txt")
findChromPeaks_milliWave
findChromPeaks_milliWave_OnDiskMSnExp
findChromPeaks_milliWave_Spectrum_list
do_findChromPeaks_milliWave
.milliWave
sink()



# Run both ----
timing <- list()
v <- Sys.time()
xdata <- findChromPeaks(raw_data, param = cwp) #Takes about 25 mins in serial on laptop
timing[[1]] <- Sys.time()-v
save(xdata, file = "xcms/xdata")
#load("xcms/xdata")

v <- Sys.time()
mdata <- findChromPeaks_milliWave(raw_data, param = cwp)
timing[[2]] <- Sys.time()-v
save(mdata, file = "xcms/milliWave_code/mdata")
#load("xcms/milliWave_code/mdata")
sink("xcms/timing.txt")
timing
sink()


# Compare structures
sink("xcms/milliWave_code/xdata_str.txt")
str(xdata)
sink()
sink("xcms/milliWave_code/mdata_str.txt")
str(mdata)
sink()

# Compare peaklists
write.csv(chromPeaks(xdata), file = "xcms/milliWave_code/xdata_peaks.csv")
write.csv(chromPeaks(mdata), file = "xcms/milliWave_code/mdata_peaks.csv")

print(dim(chromPeaks(xdata)))
print(dim(chromPeaks(mdata)))