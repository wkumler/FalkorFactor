# To demonstrate the clear superiority of Will's new noise estimation algorithms

library(xcms)
library(tidyverse)
source("xcms/SNR_improvements/will_SNRfunctions.R")

# Load relevant data ----
stds <- read.csv("xcms/Ingalls_Lab_Standards_Will.csv", stringsAsFactors = F)
good_peak <- which(stds$Compound.Name=="Leucine")

# raw_data <- "mzMLs" %>%
#   list.files(full.names = T) %>%
#   `[`(c(1,5:7,17:40)) %>%
#   readMSData(mode = "onDisk")
# save(raw_data, file = "xcms/SNR_improvements/raw_data")
load("xcms/SNR_improvements/raw_data")

# good_chr_raw <- chromatogram(raw_data,
#                              mz = c(stds$m.z[good_peak]-0.001,
#                                     stds$m.z[good_peak]+0.001),
#                              rt = c(stds$rt.sec[good_peak]-100,
#                                     stds$rt.sec[good_peak]+100))
# save(good_chr_raw, file = "xcms/SNR_improvements/good_chr_raw")
good_chr_raw <- chromatogram(raw_data,
                             mz = c(stds$m.z[good_peak]-0.001,
                                    stds$m.z[good_peak]+0.001))
#save(good_chr_raw, file = "xcms/SNR_improvements/good_chr_raw")



#load("xcms/SNR_improvements/good_chr_raw")



# Render traditional CentWave findings ----
register(SerialParam())
mock <- function(chr_raw){
  print(chromPeaks(chr_raw))
  plot(chr_raw, main="")
  par(mfrow=c(5, 6))
  par(mar=c(0.1, 0.1, 0.1, 0.1))
  for(i in chr_raw){plot(i, main="", xaxt="n", yaxt="n", ylab="", xlab="")}
  par(mfrow=c(1,1))
  par(mar=c(4.1, 4.1, 0.1, 0.1))
}

# Find peaks in all files w defaults
xchr <- findChromPeaks(good_chr_raw, param = CentWaveParam())
mock(xchr)

# Find peaks in all files with sensible defaults
xchr <- findChromPeaks(good_chr_raw, param = CentWaveParam(snthresh = 1))
mock(xchr)

# Highlight some comedy
layout(matrix(c(1,2), nrow = 2, byrow = T))
xchr <- findChromPeaks(good_chr_raw[[28]], param = CentWaveParam())
plot(xchr, main="")
abline(v=xchr@chromPeaks[1,"rt"], col="green", lwd=2)
legend("left", legend = paste("Signal-to-noise:", xchr@chromPeaks[1,"sn"]))
plot(xchr, xlim=c(xchr@chromPeaks[1, "rtmin"], xchr@chromPeaks[1, "rtmax"]),
     ylim=c(0, xchr@chromPeaks[1, "maxo"]), main="", lwd=3)
abline(v=xchr@chromPeaks[1,"rt"], col="green", lwd=3)
legend("left", legend = paste("Signal-to-noise:", xchr@chromPeaks[1,"sn"]))
layout(1)

layout(matrix(c(1,2), nrow = 2, byrow = T))
xchr <- findChromPeaks(good_chr_raw[[28]], param = CentWaveParam())
plot(xchr, main="")
abline(v=xchr@chromPeaks[3,"rt"], col="green", lwd=2)
legend("left", legend = paste("Signal-to-noise:", xchr@chromPeaks[3,"sn"]))
plot(xchr, xlim=c(xchr@chromPeaks[3, "rtmin"], xchr@chromPeaks[3, "rtmax"]),
     ylim=c(0, xchr@chromPeaks[3, "maxo"]), main="", lwd=3)
abline(v=xchr@chromPeaks[3,"rt"], col="green", lwd=3)
legend("left", legend = paste("Signal-to-noise:", xchr@chromPeaks[3,"sn"]))
layout(1)


# New algorithms
will_xchr <- will_findChromPeaks_plural(object = good_chr_raw, 
                                        param = CentWaveParam())
mock(will_xchr)
layout(matrix(c(1,2), nrow = 2, byrow = T))
will_xchr <- findChromPeaks(good_chr_raw[[28]], param = CentWaveParam())
plot(will_xchr, main="")
abline(v=will_xchr@chromPeaks[1,"rt"], col="green", lwd=2)
legend("left", legend = paste("Signal-to-noise:", will_xchr@chromPeaks[1,"sn"]))
plot(will_xchr, xlim=c(will_xchr@chromPeaks[1, "rtmin"], will_xchr@chromPeaks[1, "rtmax"]),
     ylim=c(0, will_xchr@chromPeaks[1, "maxo"]), main="", lwd=3)
abline(v=will_xchr@chromPeaks[1,"rt"], col="green", lwd=3)
legend("left", legend = paste("Signal-to-noise:", will_xchr@chromPeaks[1,"sn"]))
layout(1)
