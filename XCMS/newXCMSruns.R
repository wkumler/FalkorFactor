# XCMS, but with CAMERA this time
# Just kidding, CAMERA doesn't do what I wanted it to at ALL
# Update: xMSannotator also sucks.

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
grabSingleFileMS2 <- function(filename){
  msdata <- mzR::openMSfile(filename)
  fullhd <- mzR::header(msdata)
  ms2rows <- seq_len(nrow(fullhd))[fullhd$msLevel>1]
  spectra_list <- lapply(ms2rows, function(x){
    rtime <- fullhd[x, "retentionTime"]
    premz <- fullhd[x, "precursorMZ"]
    fragments <- mzR::peaks(msdata, x)
    return(cbind(rtime, premz, fragments))
  })
  all_data <- `names<-`(as.data.frame(do.call(rbind, spectra_list)), 
                        c("rt", "premz", "fragmz", "int"))
  return(all_data)
}
speedyQscoreCalculator <- function(file_peaks, grabSingleFileData, qscoreCalculator){
  library(data.table)
  file_data <- grabSingleFileData(unique(file_peaks$filename))
  file_data_table <- as.data.table(file_data)
  peak_splits <- split(file_peaks, ceiling(seq_len(nrow(file_peaks))/100))
  file_qscores <- lapply(peak_splits, function(i){
    eic_many <- file_data_table[mz%between%c(min(i$mzmin), max(i$mzmax))]
    individual_peaks <- split(i, seq_len(nrow(i)))
    output <- sapply(individual_peaks, function(peak_row){
      eic <- eic_many[rt%between%c(peak_row$rtmin, peak_row$rtmax) & 
                        mz%between%c(peak_row$mzmin, peak_row$mzmax)]
      qscoreCalculator(eic)
    }, USE.NAMES = FALSE)
  })
  file_peaks$sn <- unlist(file_qscores)
  return(file_peaks)
}
qscoreCalculator <- function(eic){
  #Check for bogus EICs
  if(nrow(eic)<5){
    return(0)
  }
  #Calculate where each rt would fall on a beta dist (accounts for missed scans)
  scaled_rts <- (eic$rt-min(eic$rt))/(max(eic$rt)-min(eic$rt))
  
  # Create a couple different skews and test fit
  maybe_skews <- c(2.5,3,4,5) #Add 7 to catch more multipeaks and more noise
  #Add 2 to catch very slopey peaks and more noise
  best_skew <- maybe_skews[which.max(sapply(maybe_skews, function(x){
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
  return(round(SNR*peak_cor^4*log10(max(eic$int))))
}
isIsoAdduct <- function(file_peaks, xdata, grabSingleFileData, 
                        checkPeakCor, pmppm, trapz){
  #Is the feature an isotope? I.e., is there a reasonable peak 1.003355 daltons less?
  #Load the file and apply retention time correction
  file_path <- paste("mzMLs", unique(file_peaks$file_name), sep = "/")
  message("Grabbing file data")
  file_data <- grabSingleFileData(file_path)
  message("Adjusting retention time")
  file_data$rt <- xcms::adjustedRtime(xdata)[
    MSnbase::fromFile(xdata)==unique(file_peaks$sample)][
      factor(file_data$rt)]
  library(data.table)
  file_data_table <- as.data.table(file_data)
  peak_splits <- split(file_peaks, ceiling(seq_len(nrow(file_peaks))/10))
  
  message("Beginning isotope finding")
  iso_matches_all <- lapply(peak_splits, function(i){
    eic_many <- file_data_table[mz%between%c(min(i$mzmin), max(i$mzmax))]
    individual_peaks <- split(i, seq_len(nrow(i)))
    iso_matches <- lapply(individual_peaks, function(peak_row_data){
      init_eic <- eic_many[mz%between%pmppm(peak_row_data["mz"], ppm = 5) & 
                             rt%between%c(peak_row_data["rtmin"], 
                                          peak_row_data["rtmax"])]
      init_area <- trapz(init_eic$rt, init_eic$int)
      is_M1 <- checkPeakCor(mass = peak_row_data["mz"]-1.003355,
                            rtmin=peak_row_data["rtmin"], rtmax=peak_row_data["rtmax"],
                            init_eic = init_eic, file_data_table = eic_many, 
                            pmppm = pmppm, trapz = trapz)
      is_M2 <- checkPeakCor(mass = peak_row_data["mz"]-2*1.003355,
                            rtmin=peak_row_data["rtmin"], rtmax=peak_row_data["rtmax"],
                            init_eic = init_eic, file_data_table = eic_many, 
                            pmppm = pmppm, trapz = trapz)
      is_M3 <- checkPeakCor(mass = peak_row_data["mz"]-3*1.003355,
                            rtmin=peak_row_data["rtmin"], rtmax=peak_row_data["rtmax"],
                            init_eic = init_eic, file_data_table = eic_many, 
                            pmppm = pmppm, trapz = trapz)
      is_S34 <- checkPeakCor(mass = peak_row_data["mz"]-1.995796,
                             rtmin=peak_row_data["rtmin"], rtmax=peak_row_data["rtmax"],
                             init_eic = init_eic, file_data_table = eic_many, 
                             pmppm = pmppm, trapz = trapz)
      is_Na <- checkPeakCor(mass = peak_row_data["mz"]-22.98922+1.007276,
                            rtmin=peak_row_data["rtmin"], rtmax=peak_row_data["rtmax"],
                            init_eic = init_eic, file_data_table = file_data_table, 
                            pmppm = pmppm, trapz = trapz)
      is_NH4 <- checkPeakCor(mass = peak_row_data["mz"]-18.0338+1.007276,
                             rtmin=peak_row_data["rtmin"], rtmax=peak_row_data["rtmax"],
                             init_eic = init_eic, file_data_table = file_data_table, 
                             pmppm = pmppm, trapz = trapz)
      is_H2O_H <- checkPeakCor(mass = peak_row_data["mz"]+18.0106,
                               rtmin=peak_row_data["rtmin"], rtmax=peak_row_data["rtmax"],
                               init_eic = init_eic, file_data_table = file_data_table, 
                               pmppm = pmppm, trapz = trapz)
      is_2H <- checkPeakCor(mass = peak_row_data["mz"]*2-1.007276,
                            rtmin=peak_row_data["rtmin"], rtmax=peak_row_data["rtmax"],
                            init_eic = init_eic, file_data_table = file_data_table, 
                            pmppm = pmppm, trapz = trapz)
      return(c(is_M1, is_M2, is_M3, is_S34, is_Na, is_NH4, is_H2O_H, is_2H))
    })
    return(do.call(rbind, iso_matches))
  })
  iso_matches_all <- do.call(rbind, iso_matches_all)
  colnames(iso_matches_all) <- c("M1_match", "M1_area",
                                 "M2_match", "M2_area",
                                 "M3_match", "M3_area",
                                 "S34_match", "S34_area",
                                 "Na_match", "Na_area",
                                 "NH4_match", "NH4_area",
                                 "H2O_H_match", "H2O_H_area",
                                 "2H_match", "2H_area")
  return(cbind(file_peaks, iso_matches_all))
}
checkPeakCor <- function(mass, rtmin, rtmax, init_eic, 
                         file_data_table, pmppm, trapz){
  given_eic <- file_data_table[mz%between%pmppm(mass, ppm = 5) & rt%between%c(rtmin, rtmax)]
  if(nrow(given_eic)<5){
    return(c(0, 0))
  }
  merged_eic <- merge(init_eic, given_eic, by="rt")
  if(nrow(merged_eic)<5){
    return(c(0, 0))
  }
  peak_match <- cor(merged_eic$int.x, merged_eic$int.y)
  peak_area <- trapz(merged_eic$rt, merged_eic$int.y)
  return(c(peak_match, peak_area))
}
trapz <- function(x, y) {
  m <- length(x)
  xp <- c(x, x[m:1])
  yp <- c(numeric(m), y[m:1])
  n <- 2*m
  p1 <- sum(xp[1:(n-1)]*yp[2:n]) + xp[n]*yp[1]
  p2 <- sum(xp[2:n]*yp[1:(n-1)]) + xp[1]*yp[n]
  
  return(0.5*(p1-p2))
}
# Make sure the dev version of XCMS is installed!




# Setup things ----
library(tidyverse)
library(data.table)
library(pbapply)
library(xcms)

ms_files <- "mzMLs" %>%
  list.files(pattern = ".mzML", full.names = TRUE) %>%
  normalizePath() %>%
  `[`(!grepl("Fullneg|Fullpos|QC-KM1906", x = .))
register(BPPARAM = SnowParam(tasks = length(ms_files), progressbar = TRUE))

metadata <- data.frame(
  fileid=basename(ms_files),
  sample_group=c("Blank", "Pooled", "Sample", "Std")[c(1, rep(2, 6), rep(3, 24), rep(4, 10))],
  depth=c("Blank", "Pooled", "DCM", "25m", "Std")[c(1, rep(2, 6), rep(c(rep(3, 3), rep(4, 3)), 4), rep(5, 10))],
  spindir=c("Blank", "Pooled", "Cyclone", "Anticyclone", "Std")[c(1, rep(2, 6), rep(3, 12), rep(4, 12), rep(5, 10))],
  time=c("Blank", "Pooled", "Morning", "Afternoon", "Std")[c(1, rep(2, 6), rep(c(rep(3, 6), rep(4, 6)), 2), rep(5, 10))]
) %>% new(Class = "NAnnotatedDataFrame")



# Peakpicking ----
start_time <- Sys.time()
raw_data <- readRDS(file = "XCMS/temp_data/current_raw_data.rds")
register(BPPARAM = SnowParam(tasks = length(ms_files), progressbar = TRUE))
cwp <- CentWaveParam(ppm = 2.5, peakwidth = c(15, 15), 
                     snthresh = 1, prefilter = c(0, 10000), 
                     integrate = 2, mzCenterFun = "wMean", 
                     mzdiff = 0.001, fitgauss = FALSE, 
                     noise = 5000, firstBaselineCheck = FALSE, 
                     extendLengthMSW = TRUE)
xdata <- suppressMessages(findChromPeaks(raw_data, param = cwp))
saveRDS(xdata, file = "XCMS/temp_data/current_xdata.rds")
print(Sys.time()-start_time)
# 8 minutes

start_time <- Sys.time()
xdata <- readRDS(file = "XCMS/temp_data/current_xdata.rds")
xcms_peakdf <- chromPeaks(xdata) %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  mutate(filename=fileNames(xdata)[.[["sample"]]]) %>%
  arrange(mz)
split_xcms_filedata <- split(xcms_peakdf, xcms_peakdf$filename)
files_qscores <- bplapply(X = split_xcms_filedata, FUN = speedyQscoreCalculator,
                          grabSingleFileData, qscoreCalculator)

threshold <- 25
peakdf_qscored <- files_qscores %>%
  do.call(what = rbind) %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  filter(sn>threshold) %>%
  select(-filename) %>%
  arrange(sample, rtmin, rtmax) %>%
  as.matrix()
write.csv(peakdf_qscored, file = "XCMS/temp_data/peakdf_qscored.csv", 
          row.names = FALSE)
xdata_cleanpeak <- `chromPeaks<-`(xdata, peakdf_qscored)
saveRDS(xdata_cleanpeak, file = "XCMS/temp_data/xdata_cleanpeak.rds")
print(Sys.time()-start_time)
# 10 minutes?



# Plot picked peaks, just for kicks! ----
xdata_cleanpeak <- readRDS("XCMS/temp_data/xdata_cleanpeak.rds")
given_file <- ms_files[2]
clean_peakdf <- chromPeaks(xdata_cleanpeak) %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  mutate(filename=fileNames(xdata_cleanpeak)[.[["sample"]]]) %>%
  arrange(mz) %>%
  filter(filename==given_file)
file_eic <- as.data.table(grabSingleFileData(given_file))
pdf(width = 17, height = 11, file = "XCMS")
par(mfrow=c(6, 8))
par(mar=c(2.1, 2.1, 0.1, 0.1))
apply(clean_peakdf, 1, function(x){
  eic <- file_eic[rt%between%(as.numeric(c(x["rtmin"], x["rtmax"]))+c(-50, 50))&
                    mz%between%as.numeric(c(x["mzmin"], x["mzmax"]))]
  plot(eic$rt, eic$int, type="l", xlab="", ylab="", 
       ylim=c(0, max(eic$int)*1.2), xlim=c(min(eic$rt)-10, max(eic$rt)+10))
  abline(v=as.numeric(x["rtmin"]), col="red")
  abline(v=as.numeric(x["rtmax"]), col="red")
  legend("topright", legend = x["sn"], bty="n")
  legend("topleft", legend = substr(x["mz"], 1, 9), bty="n")
  return(NULL)
})
dev.off()



# Other XCMS things (rtcor, group) ----
start_time <- Sys.time()
register(BPPARAM = SerialParam())
register(BPPARAM = SnowParam(tasks = length(ms_files), progressbar = TRUE))
obp <- ObiwarpParam(binSize = 0.1, centerSample = 4, 
                    response = 1, distFun = "cor_opt")
xdata_rt <- suppressMessages(adjustRtime(xdata_cleanpeak, param = obp))
plotAdjustedRtime(xdata_rt, col = c("green", "red", "blue", "black", "black")[
  factor(metadata@data$depth)])

pdp <- PeakDensityParam(sampleGroups = xdata_rt$depth, 
                        bw = 5, minFraction = 0.5, 
                        binSize = 0.002, minSamples = 2)
xdata_cor <- groupChromPeaks(xdata_rt, param = pdp)

fpp <- FillChromPeaksParam()
xdata_filled <- suppressMessages(fillChromPeaks(xdata_cor, param = fpp))

saveRDS(xdata_filled, file = "XCMS/temp_data/current_xdata_filled.rds")
print(Sys.time()-start_time)
# 10 minutes

show(featureDefinitions(xdata_filled))



# Find isotopes and adducts ----
xdata_filled <- readRDS("XCMS/temp_data/current_xdata_filled.rds")
feature_defs <- featureDefinitions(xdata_filled)
feature_peaks <- lapply(seq_len(nrow(feature_defs)), function(i){
  cbind(feature=sprintf("FT%03d", i), 
        peak_id=unlist(feature_defs$peakidx[i]))
}) %>% 
  do.call(what=rbind) %>% 
  as.data.frame(stringsAsFactors=FALSE) %>% 
  mutate(peak_id=as.numeric(peak_id)) %>%
  cbind(chromPeaks(xdata_filled)[.$peak_id, ]) %>%
  mutate(file_name=basename(fileNames(xdata_filled))[sample]) %>%
  arrange(feature, sample)

is_peak_iso <- bplapply(split(feature_peaks, feature_peaks$file_name), 
                        FUN = isIsoAdduct, xdata=xdata_filled,
                        grabSingleFileData=grabSingleFileData,
                        checkPeakCor=checkPeakCor, 
                        pmppm=pmppm, trapz=trapz) %>%
  do.call(what = rbind) %>% as.data.frame()

peakshapematch_features <- is_peak_iso %>%
  group_by(feature) %>%
  summarise(prob_M1=median(M1_match), prob_M2=median(M2_match), prob_M3=median(M3_match),
            prob_Na=median(Na_match), prob_NH4=median(NH4_match), 
            prob_H2O_H=median(H2O_H_match), prob_2H=median(`2H_match`)) %>%
  rowwise() %>%
  mutate(is_funny=max(prob_M1, prob_M2, prob_M3)) %>%
  ungroup() %>%
  filter(is_funny>0.8)

responsematch_features <- is_peak_iso %>%
  

possible_M1s <- goodcor_features %>%
  filter(prob_M1>0.8) %>%
  pull(feature) %>%
  `[`(xcms::featureDefinitions(xdata_filled), .,) %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  select(mzmin, mzmax, rtmin, rtmax) %>%
  split(seq_len(nrow(.)))



removed_features[,-1] <- round(removed_features[,-1], digits = 4)
removed_features <- as.data.frame(removed_features)
removed_features[removed_features<0.8] <- "-------"
removed_features
