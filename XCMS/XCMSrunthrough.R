
# Setup things ----
library(xcms)
library(dplyr)
library(RSQLite)
library(data.table)
library(beepr)
library(pbapply)
# if(dir.exists("temp_data"))unlink("temp_data")
# if(!dir.exists("temp_data"))dir.create("temp_data")
start_time <- Sys.time()



# Functions ----
pmppm <- function(mass, ppm=4){c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))}
grabSingleFileData <- function(filename){
  msdata <- mzR:::openMSfile(filename)
  fullhd <- mzR::header(msdata)
  spectra_list <- lapply(seq_len(nrow(fullhd)), function(x){
    given_peaks <- mzR::peaks(msdata, x)
    rtime <- fullhd[x, "retentionTime"]
    return(cbind(rtime, given_peaks))
  })
  all_data <- `names<-`(as.data.frame(do.call(rbind, spectra_list)), 
                        c("rt", "mz", "int"))
  return(all_data)
}
qscoreCalculator <- function(eic){
  #Check for bogus EICs
  if(nrow(eic)<5){
    return(0)
  }
  #Create an "ideal" peak of the same width
  perf_peak <- dnorm(seq(-3, 3, length.out = nrow(eic)))
  #Calculate the correlation between the perfect peak and the observed values
  peak_cor <- cor(perf_peak, eic$int)
  #Calculate the normalized residuals
  residuals <- eic$int/max(eic$int)-perf_peak/max(perf_peak)
  #Calculate the minimum SD, after normalizing for any shape discrepancy
  old_res_sd <- sd(residuals)
  norm_residuals <- diff(residuals)
  new_res_sd <- sd(norm_residuals)
  while(new_res_sd<old_res_sd){
    old_res_sd <- new_res_sd
    norm_residuals <- diff(residuals)
    new_res_sd <- sd(residuals)
  }
  #Calculate SNR
  SNR <- (max(eic$int)-min(eic$int))/sd(norm_residuals*max(eic$int))
  #Return the quality score
  output <- data.frame(SNR, peak_cor, qscore=SNR*peak_cor^4*log10(max(eic$int)))
  return(output)
}
xcmsQscoreCalculator <- function(df_row, xcms_peakdf, file_data_table){
  #Extract the relevant EIC
  peak_row_data <- xcms_peakdf[df_row, ]
  eic <- file_data_table[rt %between% c(peak_row_data$rtmin, peak_row_data$rtmax)&
                           mz %between% c(peak_row_data$mzmin, peak_row_data$mzmax)]
  return(qscoreCalculator(eic))
}

# Load MS data ----
ms_files <- "mzMLs" %>%
  list.files(pattern = ".mzML", full.names = TRUE) %>%
  normalizePath() %>%
  `[`(!grepl("Fullneg|Fullpos|QC-KM1906", x = .))

metadata <- data.frame(
  fileid=1:41, 
  filenames=gsub("190715_", "", basename(ms_files)),
  sample_group=c("Blank", "Pooled", "Sample", "Std")[c(1, rep(2, 6), rep(3, 24), rep(4, 10))],
  depth=c("Blank", "Pooled", "DCM", "25m", "Std")[c(1, rep(2, 6), rep(c(rep(3, 3), rep(4, 3)), 4), rep(5, 10))],
  spindir=c("Blank", "Pooled", "Cyclone", "Anticyclone", "Std")[c(1, rep(2, 6), rep(3, 12), rep(4, 12), rep(5, 10))],
  time=c("Blank", "Pooled", "Morning", "Afternoon", "Std")[c(1, rep(2, 6), rep(c(rep(3, 6), rep(4, 6)), 2), rep(5, 10))]
) %>%
  new(Class = "NAnnotatedDataFrame")

raw_data <- readMSData(files = ms_files, pdata = metadata, mode = "onDisk")
saveRDS(raw_data, file = "temp_data/current_raw_data.rds")
print(Sys.time()-start_time)
#2.23 minutes
beep(2)



# Perform peakpicking ----
raw_data <- readRDS(file = "XCMS/temp_data/current_raw_data.rds")
cwp <- CentWaveParam(ppm = 5, peakwidth = c(20, 80), 
                     snthresh = 0, prefilter = c(0, 0), 
                     integrate = 1, mzCenterFun = "wMean", 
                     mzdiff = 0.0001, fitgauss = FALSE, 
                     noise = 0, firstBaselineCheck = FALSE)
xdata <- findChromPeaks(raw_data, param = cwp)
print(xdata)
saveRDS(xdata, file = "XCMS/temp_data/current_xdata.rds")
print(Sys.time()-start_time)
# 34.2 minutes
beep(2)



# Re-assign quality scores to confirm good peaks ----
xdata <- readRDS(file = "XCMS/temp_data/current_xdata.rds")
xcms_peakdf <- as.data.frame(chromPeaks(xdata))

fileids <- unique(xcms_peakdf$sample)
split_xcms_filedata <- split(xcms_peakdf, xcms_peakdf$sample)
start_time <- Sys.time()
files_qscores <- bplapply(fileids, function(x){
  print(paste("Processing", ms_files[x]))
  file_peaks <- split_xcms_filedata[[x]]
  file_data <- grabSingleFileData(ms_files[x])
  file_data_table <- as.data.table(file_data)
  file_qscores <- lapply(seq_len(nrow(file_peaks)), 
                           FUN = xcmsQscoreCalculator, 
                           xcms_peakdf=file_peaks, 
                           file_data_table=file_data_table)
  file_qscores_df <- do.call(rbind, file_qscores)
  return(cbind(split_xcms_filedata[[x]], file_qscores_df))
}, SnowParam(progressbar = TRUE, tasks = length(fileids)))
peakdf_qscored <- as.data.frame(do.call(rbind, files_qscores))
write.csv(peakdf_qscored, file = "XCMS/temp_data/peakdf_qscored.csv", row.names = FALSE)
print(Sys.time()-start_time)
beep(2)



# Decide on quality threshold and reassign peaklist ----
peakdf_qscored <- read.csv("XCMS/temp_data/peakdf_qscored.csv")
threshold <- 20
cleandf_qscored <- peakdf_qscored[peakdf_qscored$qscore>threshold,]
xdata_cleanpeak <- `chromPeaks<-`(xdata, cleandf_qscored)



# Adjust retention time and compare ----
obp <- ObiwarpParam(binSize = 0.01, centerSample = 4, response = 1, distFun = "cor_opt")
xdata_rt <- adjustRtime(xdata_cleanpeak, param = obp)
plotAdjustedRtime(xdata_rt)
saveRDS(xdata_rt, file = "temp_data/current_xdata_rt.rds")
print(Sys.time()-start_time)
# 50 minutes
beep(2)



# Correspondence ----
xdata_rt <- readRDS("temp_data/current_xdata_rt.rds")
pdp <- PeakDensityParam(sampleGroups = xdata_rt$sample_group, 
                        bw = 20, minFraction = 0.5, 
                        binSize = 0.002)
xdata_cor <- groupChromPeaks(xdata_rt, param = pdp)
print(Sys.time()-start_time)
beep(2)


# Group isotopes and adducts ----

