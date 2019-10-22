#microWave heart demo!

library(xcms)
library(tidyverse)
library(roxygen2)

#source("microWave_functions.R")
source("xcms/microWave_code/microWave_functions.R")

# msfiles <- list.files(path = "mzMLs", pattern = "Blk|Smp|Full\\d", full.names = T)
# raw_data <- readMSData(files = msfiles, msLevel. = 1, centroided. = T, mode = "onDisk")
# save(raw_data, file = "xcms/raw_data")
load("xcms/raw_data")
x <- filterMsLevel(raw_data, msLevel. = 1L)
x <- selectFeatureData(x, fcol = c(MSnbase:::.MSnExpReqFvarLabels, "centroided"))
x <- lapply(1:length(fileNames(x)), FUN=filterFile, object = x)
x <- x[[2]]
x <- spectra(x)

all_data <- data.frame(mz=unlist(lapply(x, mz), use.names = FALSE),
                       int=unlist(lapply(x, intensity), use.names = FALSE), 
                       rt=rep(unname(unlist(lapply(x, rtime))), sapply(lapply(x, mz), length)))


data <- all_data %>% filter(mz>100&mz<120) %>% filter(rt>60&rt<1100)

eic_list <- constructEICs(data)

peak_df <- microWavePeaks(eic_list)
peak_df <- arrange(peak_df, desc(Peak_SNR*Peak_gauss_fit^4))

peakCheck(eic_list, peak_df, "1.1.1")
for(i in peak_df$Peak_id){
  peakCheck(eic_list, peak_df, i)
  replot <- readline(prompt = "Press Enter") 
  if(replot=="j"){
    peakCheck(eic_list, peak_df, i, zoom=T)
    readline(prompt = "Continue?")
  } else if(replot=="k") {
    peakCheck(eic_list, peak_df, i, pts = T, zoom = T)
    readline(prompt = "Continue?")
  } else {
    next
  }
}
