# Code for collecting source files of various CentWave algorithms
load("xcms/chr_raw")

# Applied to chr_raw, a Chromatograms object
showMethods("findChromPeaks")
findMethods("findChromPeaks")[[3]]
#Basically just a wrapper for .findChromPeaks_XChromatograms()
#Activates bpparam()


# Sourcing .findChromPeaks_XChromatograms()
xcms:::.findChromPeaks_XChromatograms
#Applies findChromPeaks in parallel via bplapply to each chromatogram in chr_raw
#lapply(chr_raw@.Data, print)


# Applied to a single Chromatogram in chr_raw
showMethods("findChromPeaks")
findMethods("findChromPeaks")[[1]]
#Applies peaksWithCentWave to various sections of chr_raw
# I.e. accepts int, rt, and param
#Fills in the chromPeaks() slot of the Chromatogram object


# Applying peaksWithCentWave to a single Chromatogram
int <- intensity(chr_raw@.Data[[1]])
rt <- rtime(chr_raw@.Data[[1]])
plot(int~rt)
peak1data <- peaksWithCentWave(int, rt)
abline(v=peak1data[colnames(peak1data)=="rtmin"])
abline(v=peak1data[colnames(peak1data)=="rtmax"])
abline(v=peak1data[colnames(peak1data)=="intb"])


# Looking into peaksWithCentwave when applied to single Chromatogram
sink("will_peaksWithCentWave.R")
peaksWithCentWave
sink()
readLines("will_peaksWithCentWave.R") %>% 
  c("will_peaksWithCentWave <- ", .) %>%
  .[-c(length(.)-1, length(.))] %>%
  gsub(pattern = ".getRtROI", replacement = "xcms:::.getRtROI") %>%
  gsub(pattern = "continuousPtsAboveThreshold", replacement = "xcms:::continuousPtsAboveThreshold") %>%
  gsub(pattern = "joinOverlappingPeaks", replacement = "xcms:::joinOverlappingPeaks") %>%
  gsub(pattern = ".narrow_rt_boundaries", replacement = "xcms:::.narrow_rt_boundaries") %>%
  gsub(pattern = "rectUnique", replacement = "xcms:::rectUnique") %>%
  writeLines("will_peaksWithCentWave.R")
source("will_peaksWithCentWave.R")
peak1data <- will_peaksWithCentWave(int, rt)
abline(v=peak1data[colnames(peak1data)=="rtmin"])
abline(v=peak1data[colnames(peak1data)=="rtmax"])
abline(v=peak1data[colnames(peak1data)=="intb"])

