#PeakCheck
library(tidyverse)
library(xcms)
library(data.table)
library(plotly)

# Functions ----
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
    return(data.frame(SNR=0, peak_cor=0, qscore=0))
  }
  #Calculate where each rt would fall on a beta distribution (accounts for missed scans)
  scaled_rts <- (eic$rt-min(eic$rt))/(max(eic$rt)-min(eic$rt))
  
  # Create a couple different skews and test fit
  possible_skews <- c(2,3,5)
  best_skew <- possible_skews[which.max(sapply(possible_skews, function(x){
    cor(dbeta(scaled_rts, shape1 = x, shape2 = 5), eic$int)
  }))]
  perf_peak <- dbeta(scaled_rts, shape1 = best_skew, shape2 = 5)
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
xcmsQscoreCalculator <- function(df_row, xcms_peakdf, file_data_table, 
                                 qscoreCalculator = qscoreCalculator){
  #Extract the relevant EIC
  peak_row_data <- xcms_peakdf[df_row, ]
  print(peak_row_data)
  eic <- file_data_table[rt %between% c(peak_row_data$rtmin, peak_row_data$rtmax)&
                           mz %between% c(peak_row_data$mzmin, peak_row_data$mzmax)]
  return(qscoreCalculator(eic))
}



# Other ----

xdata <- readRDS(file = "XCMS/temp_data/current_xdata.rds")
xdata_filled <- readRDS(file = "XCMS/temp_data/current_xdata_filled.rds")
peakdf_qscored <- as.matrix(read.csv(file = "XCMS/temp_data/peakdf_qscored.csv"))
ms_files <- "mzMLs" %>%
  list.files(pattern = ".mzML", full.names = TRUE) %>%
  normalizePath() %>%
  `[`(!grepl("Fullneg|Fullpos|QC-KM1906", x = .))
metadata <- data.frame(
  fileid=seq_along(ms_files), 
  filenames=gsub("190715_", "", basename(ms_files)),
  sample_group=c("Blank", "Pooled", "Sample", "Std")[c(1, rep(2, 6), rep(3, 24), rep(4, 10))],
  depth=c("Blank", "Pooled", "DCM", "25m", "Std")[c(1, rep(2, 6), rep(c(rep(3, 3), rep(4, 3)), 4), rep(5, 10))],
  spindir=c("Blank", "Pooled", "Cyclone", "Anticyclone", "Std")[c(1, rep(2, 6), rep(3, 12), rep(4, 12), rep(5, 10))],
  time=c("Blank", "Pooled", "Morning", "Afternoon", "Std")[c(1, rep(2, 6), rep(c(rep(3, 6), rep(4, 6)), 2), rep(5, 10))]
)



all_file_MS1 <- pblapply(seq_along(ms_files), function(x){
  file_data <- grabSingleFileData(ms_files[x])
  file_data$rt <- xcms::adjustedRtime(xdata_filled)[
    MSnbase::fromFile(xdata_filled)==x][factor(file_data$rt)]
  as.data.table(cbind(fileid=x, file_data))
}) %>% do.call(what = rbind)

missed_peak <- 133.037509+1.007276 #Aspartic acid
missed_peak <- 175.095691+1.007276 #Citrulline
eic <- all_file_MS1[mz%between%pmppm(missed_peak, ppm = 5)]
eic_metadata <- left_join(eic, metadata, by="fileid")
gp <- ggplot(eic_metadata) + geom_line(aes(x=rt, y=int, group=fileid, 
                                           color=depth, text=depth)) +
  scale_color_manual(values = c("25m"="#0000FF99", "Blank"="#FF000099", 
                                "DCM"="#00FF0099", "Pooled"="#00000099",
                                "Std"="#00000099"), 
                     name="Depth") + theme_bw()
ggplotly(gp)

given_peakdf <- peakdf_qscored %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  filter(mz%between%pmppm(missed_peak, ppm = 5))

for(i in seq_len(nrow(given_peakdf))){
  peak_row_data <- given_peakdf[i,]
  eic <- all_file_MS1[rt %between% c(peak_row_data$rtmin, peak_row_data$rtmax)&
                        mz %between% c(peak_row_data$mzmin, peak_row_data$mzmax)&
                        fileid==i]
  print(qscoreCalculator(eic))
}


par(mfrow=c(3,3))
