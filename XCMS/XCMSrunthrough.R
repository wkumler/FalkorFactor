
### Setup things ----
library(xcms)
library(tidyverse)
library(data.table)
library(beepr)
library(httr)
library(pbapply)
library(Rdisop)
library(fastmatch)
start_time <- Sys.time()

ms_files <- "mzMLs" %>%
  list.files(pattern = ".mzML", full.names = TRUE) %>%
  normalizePath() %>%
  `[`(!grepl("Fullneg|Fullpos|QC-KM1906", x = .))
register(BPPARAM = SnowParam(tasks = length(ms_files), progressbar = TRUE))

metadata <- data.frame(
  fileid=seq_along(ms_files), 
  filenames=gsub("190715_", "", basename(ms_files)),
  sample_group=c("Blank", "Pooled", "Sample", "Std")[c(1, rep(2, 6), rep(3, 24), rep(4, 10))],
  depth=c("Blank", "Pooled", "DCM", "25m", "Std")[c(1, rep(2, 6), rep(c(rep(3, 3), rep(4, 3)), 4), rep(5, 10))],
  spindir=c("Blank", "Pooled", "Cyclone", "Anticyclone", "Std")[c(1, rep(2, 6), rep(3, 12), rep(4, 12), rep(5, 10))],
  time=c("Blank", "Pooled", "Morning", "Afternoon", "Std")[c(1, rep(2, 6), rep(c(rep(3, 6), rep(4, 6)), 2), rep(5, 10))]
) %>% new(Class = "NAnnotatedDataFrame")



### Functions ----
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
  eic <- file_data_table[rt %between% c(peak_row_data$rtmin, peak_row_data$rtmax)&
                           mz %between% c(peak_row_data$mzmin, peak_row_data$mzmax)]
  return(qscoreCalculator(eic))
}
isIso <- function(file_peaks, xdata, grabSingleFileData, checkPeakCor, pmppm){
  #Is the feature an isotope? I.e., is there a reasonable peak 1.003355 daltons less?
  #Load the file and apply retention time correction
  file_path <- paste("mzMLs", unique(file_peaks$file_name), sep = "/")
  file_data <- grabSingleFileData(file_path)
  file_data$rt <- xcms::adjustedRtime(xdata)[
    MSnbase::fromFile(xdata)==unique(file_peaks$sample)][
      factor(file_data$rt)]
  library(data.table)
  file_dt <- as.data.table(file_data)
  
  iso_matches <- t(apply(file_peaks[,c("mz","rtmin","rtmax")], 1, function(peak_row_data){
    init_eic <- file_dt[mz%between%pmppm(peak_row_data["mz"], ppm = 5) & 
                          rt%between%c(peak_row_data["rtmin"], 
                                       peak_row_data["rtmax"])]
    is_M1 <- checkPeakCor(mass = peak_row_data["mz"]-1.003355,
                          rtmin=peak_row_data["rtmin"], rtmax=peak_row_data["rtmax"],
                          init_eic = init_eic, file_dt = file_dt, pmppm = pmppm)
    is_M2 <- checkPeakCor(mass = peak_row_data["mz"]-2*1.003355,
                          rtmin=peak_row_data["rtmin"], rtmax=peak_row_data["rtmax"],
                          init_eic = init_eic, file_dt = file_dt, pmppm = pmppm)
    is_M3 <- checkPeakCor(mass = peak_row_data["mz"]-3*1.003355,
                          rtmin=peak_row_data["rtmin"], rtmax=peak_row_data["rtmax"],
                          init_eic = init_eic, file_dt = file_dt, pmppm = pmppm)
    is_S34 <- checkPeakCor(mass = peak_row_data["mz"]-1.995796,
                          rtmin=peak_row_data["rtmin"], rtmax=peak_row_data["rtmax"],
                          init_eic = init_eic, file_dt = file_dt, pmppm = pmppm)
    
    return(cbind(is_M1, is_M2, is_M3, is_S34))
  }))
  colnames(iso_matches) <- c("M1_match", "M2_match", "M3_match", "S34_match")
  return(cbind(file_peaks, iso_matches))
}
isAdduct <- function(file_peaks, xdata, grabSingleFileData, checkPeakCor, pmppm){
  #Is the feature an adduct? I.e., is there a reasonable peak at the [M+H] mass too?
  file_path <- paste("mzMLs", unique(file_peaks$file_name), sep = "/")
  file_data <- grabSingleFileData(file_path)
  file_data$rt <- xcms::adjustedRtime(xdata)[
    MSnbase::fromFile(xdata)==unique(file_peaks$sample)][
      factor(file_data$rt)]
  library(data.table)
  file_dt <- as.data.table(file_data)
  
  adduct_matches <- t(apply(file_peaks[,c("mz","rtmin","rtmax")], 1, function(peak_row_data){
    init_eic <- file_dt[mz%between%pmppm(peak_row_data["mz"], ppm = 5) & 
                          rt%between%c(peak_row_data["rtmin"], 
                                       peak_row_data["rtmax"])]
    is_Na <- checkPeakCor(mass = peak_row_data["mz"]-22.98922+1.007276,
                          rtmin=peak_row_data["rtmin"], rtmax=peak_row_data["rtmax"],
                          init_eic = init_eic, file_dt = file_dt, pmppm = pmppm)
    is_NH4 <- checkPeakCor(mass = peak_row_data["mz"]-18.0338+1.007276,
                          rtmin=peak_row_data["rtmin"], rtmax=peak_row_data["rtmax"],
                          init_eic = init_eic, file_dt = file_dt, pmppm = pmppm)
    is_H2O_H <- checkPeakCor(mass = peak_row_data["mz"]+18.0106,
                             rtmin=peak_row_data["rtmin"], rtmax=peak_row_data["rtmax"],
                             init_eic = init_eic, file_dt = file_dt, pmppm = pmppm)
    is_2H <- checkPeakCor(mass = peak_row_data["mz"]*2-1.007276,
                          rtmin=peak_row_data["rtmin"], rtmax=peak_row_data["rtmax"],
                          init_eic = init_eic, file_dt = file_dt, pmppm = pmppm)
    
    return(cbind(is_Na, is_NH4, is_H2O_H, is_2H))
  }))
  colnames(adduct_matches) <- c("Na_match", "NH4_match", "H2O_H_match", "2H_match")
  return(cbind(file_peaks, adduct_matches))
}
checkPeakCor <- function(mass, rtmin, rtmax, init_eic, file_dt, pmppm){
  given_eic <- file_dt[mz%between%pmppm(mass, ppm = 5) & rt%between%c(rtmin, rtmax)]
  if(nrow(given_eic)<5){
    return(0)
  }
  merged_eic <- merge(init_eic, given_eic, by="rt")
  if(nrow(merged_eic)<5){
    return(0)
  }
  peak_match <- cor(merged_eic$int.x, merged_eic$int.y)
  return(peak_match)
}
findIsoAdds <- function(file_peaks, xdata, grabSingleFileData, 
                        getPeakArea, pmppm, trapz){
  file_path <- paste("mzMLs", unique(file_peaks$file_name), sep = "/")
  file_data <- grabSingleFileData(file_path)
  file_data$rt <- xcms::adjustedRtime(xdata)[
    MSnbase::fromFile(xdata)==unique(file_peaks$sample)][
      factor(file_data$rt)]
  library(data.table)
  file_dt <- as.data.table(file_data)
  
  adduct_matches <- t(apply(file_peaks[,c("mz","rtmin","rtmax")], 1, function(peak_row_data){
    init_eic <- file_dt[mz%between%pmppm(peak_row_data["mz"], ppm = 5) & 
                          rt%between%c(peak_row_data["rtmin"], 
                                       peak_row_data["rtmax"])]
    isoadd_masses <- as.numeric(peak_row_data["mz"])+
      c(M_1=1.003355, M_2=2*1.003355, M_3=3*1.003355, S34=1.995796,
        Na=22.98922-1.007276, NH4=18.0338-1.007276,
        H2O_H=-18.0106, X2H=as.numeric(peak_row_data["mz"])-1.007276)
    isoadd_matches <- lapply(isoadd_masses, getPeakArea, 
                             rtmin=peak_row_data["rtmin"], rtmax=peak_row_data["rtmax"],
                             init_eic = init_eic, file_dt = file_dt, 
                             pmppm = pmppm, trapz = trapz)
    return(unlist(isoadd_matches))
  }))
  colnames(adduct_matches) <- c("M_1_match", "M_1_area", 
                                "M_2_match", "M_2_area",
                                "M_3_match", "M_3_area",
                                "S_34_match", "S_34_area",
                                "Na_match", "Na_area",
                                "NH4_match", "NH4_area", 
                                "H2O_H_match", "H2O_H_area",
                                "X2H_match", "X2H_area")
  return(cbind(file_peaks, adduct_matches))
}
getPeakArea <- function(mass, rtmin, rtmax, init_eic, file_dt, pmppm, trapz){
  given_eic <- file_dt[mz%between%pmppm(mass, ppm = 5) & rt%between%c(rtmin, rtmax)]
  if(nrow(given_eic)<5){
    return(c(match=0, area=0))
  }
  merged_eic <- merge(init_eic, given_eic, by="rt")
  if(nrow(merged_eic)<5){
    return(c(match=0, area=0))
  }
  peak_match <- cor(merged_eic$int.x, merged_eic$int.y)
  peak_area <- trapz(merged_eic$rt, merged_eic$int.y)
  return(c(match=peak_match, area=peak_area))
}
goldenRules <- function(formula){
  elem_data <- "[[:upper:]][[:lower:]]*[[:digit:]]*" %>%
    gregexpr(text = formula) %>%
    regmatches(x = formula) %>%
    unlist()
  single_elems <- regexpr("[[:digit:]]$", elem_data) == -1
  elem_data[single_elems] <- paste0(elem_data[single_elems], 1)
  
  elem_names <- regmatches(elem_data, regexpr("[[:alpha:]]+", elem_data))
  elem_counts <- regmatches(elem_data, regexpr("[[:digit:]]+", elem_data)) %>%
    as.numeric() %>% `names<-`(elem_names)
  #Basic SENIOR check (Rule 2, S_i only)
  senior_i <- sum(elem_counts[c("H", "N", "P")], na.rm = TRUE)%%2==0
  
  #Element ratios (Rules 4 and 5)
  mins <- elem_counts[c("H", "O", "N", "P", "S")]/elem_counts["C"]>c(0.2, numeric(4))
  maxs <- elem_counts[c("H", "O", "N", "P", "S")]/elem_counts["C"]<c(3.1, 1.2, 1.3, 0.3, 0.8)
  rule_4_check <- all(mins, maxs, na.rm = TRUE)
  
  #Element balance (Rule 6)
  # if(all(elem_counts[c("N", "O", "P", "S")])>1){
  #   rule_6_check <- all(elem_counts[c("N", "O", "P", "S")]<c(10, 20, 4, 3), na.rm=TRUE)
  # } #Nope, idk how to program this one.
  
  return(all(senior_i, rule_4_check))
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


### Load MS data ----
ms_files <- "mzMLs" %>%
  list.files(pattern = ".mzML", full.names = TRUE) %>%
  normalizePath() %>%
  `[`(!grepl("Fullneg|Fullpos|QC-KM1906", x = .))

raw_data <- readMSData(files = ms_files, pdata = metadata, 
                       mode = "onDisk", verbose = TRUE)
saveRDS(raw_data, file = "XCMS/temp_data/current_raw_data.rds")
print(Sys.time()-start_time)
#2.3 minutes
beep(2)



### Perform peakpicking ----
start_time <- Sys.time()
raw_data <- readRDS(file = "XCMS/temp_data/current_raw_data.rds")
register(BPPARAM = SnowParam(tasks = length(ms_files), progressbar = TRUE))
cwp <- CentWaveParam(ppm = 5, peakwidth = c(20, 80), 
                     snthresh = 0, prefilter = c(0, 0), 
                     integrate = 1, mzCenterFun = "wMean", 
                     mzdiff = 0.0001, fitgauss = FALSE, 
                     noise = 0, firstBaselineCheck = FALSE)
xdata <- suppressMessages(findChromPeaks(raw_data, param = cwp))
saveRDS(xdata, file = "XCMS/temp_data/current_xdata.rds")
print(Sys.time()-start_time)
# 30 minutes
beep(2)



### Re-assign quality scores to confirm good peaks ----
start_time <- Sys.time()
xdata <- readRDS(file = "XCMS/temp_data/current_xdata.rds")
xcms_peakdf <- as.data.frame(chromPeaks(xdata))

fileids <- unique(xcms_peakdf$sample)
split_xcms_filedata <- split(xcms_peakdf, xcms_peakdf$sample)
files_qscores <- bplapply(fileids, function(x, split_xcms_filedata, ms_files, 
                                            grabSingleFileData, xcmsQscoreCalculator,
                                            qscoreCalculator){
  library(data.table)
  file_peaks <- split_xcms_filedata[[x]]
  file_data <- grabSingleFileData(ms_files[x])
  file_data_table <- as.data.table(file_data)
  file_qscores <- lapply(seq_len(nrow(file_peaks)), FUN = xcmsQscoreCalculator, 
                         xcms_peakdf=file_peaks, file_data_table=file_data_table,
                         qscoreCalculator=qscoreCalculator)
  file_qscores_df <- do.call(rbind, file_qscores)
  return(cbind(split_xcms_filedata[[x]], file_qscores_df))
}, 
split_xcms_filedata=split_xcms_filedata, 
ms_files=ms_files, 
grabSingleFileData = grabSingleFileData,
xcmsQscoreCalculator = xcmsQscoreCalculator,
qscoreCalculator = qscoreCalculator)
peakdf_qscored <- as.data.frame(do.call(rbind, files_qscores))
write.csv(peakdf_qscored, file = "XCMS/temp_data/peakdf_qscored.csv", row.names = FALSE)
print(Sys.time()-start_time)
# 1.5 hours
beep(2)



### Decide on quality threshold and reassign peaklist ----
start_time <- Sys.time()
xdata <- readRDS(file = "XCMS/temp_data/current_xdata.rds")
peakdf_qscored <- as.matrix(read.csv(file = "XCMS/temp_data/peakdf_qscored.csv"))
threshold <- 50
cleandf_qscored <- peakdf_qscored[peakdf_qscored[, "qscore"]>threshold,]
xdata_cleanpeak <- `chromPeaks<-`(xdata, cleandf_qscored)
# sample_df <- peakdf_qscored[sample(1:nrow(peakdf_qscored), size = 10000),]
# plot(log10(sample_df[,"sn"]), log10(sample_df[,"qscore"]))



### Adjust retention time and compare ----
start_time <- Sys.time()
obp <- ObiwarpParam(binSize = 0.01, centerSample = 4, response = 1, distFun = "cor_opt")
xdata_rt <- suppressMessages(adjustRtime(xdata_cleanpeak, param = obp))
plotAdjustedRtime(xdata_rt, col = c("green", "red", "blue", "black", "black")[
  factor(factor(metadata@data$depth))
])
saveRDS(xdata_rt, file = "XCMS/temp_data/current_xdata_rt.rds")
print(Sys.time()-start_time)
# 17 minutes
beep(2)



### Correspondence ----
start_time <- Sys.time()
xdata_rt <- readRDS(file = "XCMS/temp_data/current_xdata_rt.rds")
pdp <- PeakDensityParam(sampleGroups = xdata_rt$depth, 
                        bw = 2, minFraction = 0.5, 
                        binSize = 0.002)
xdata_cor <- groupChromPeaks(xdata_rt, param = pdp)
saveRDS(xdata_cor, file = "XCMS/temp_data/current_xdata_cor.rds")
print(Sys.time()-start_time)
# 5 seconds
beep(2)



### Fill peaks ----
start_time <- Sys.time()
xdata_cor <- readRDS(file = "XCMS/temp_data/current_xdata_cor.rds")
xdata_filled <- suppressMessages(fillChromPeaks(xdata_cor, param = FillChromPeaksParam()))
saveRDS(object = xdata_filled, file = "XCMS/temp_data/current_xdata_filled.rds")
print(Sys.time()-start_time)
# 1.5 minutes
beep(2)


### Add new quality scores and set retention time limits ----
start_time <- Sys.time()
xdata_filled <- readRDS(file = "XCMS/temp_data/current_xdata_filled.rds")
feature_peaks <- lapply(seq_len(nrow(featureDefinitions(xdata_filled))), function(i){
  cbind(feature=sprintf("FT%03d", i), 
        peak_id=unlist(featureDefinitions(xdata_filled)$peakidx[i]))
}) %>% 
  do.call(what=rbind) %>% 
  as.data.frame(stringsAsFactors=FALSE) %>% 
  mutate(peak_id=as.numeric(peak_id)) %>%
  cbind(chromPeaks(xdata_filled)[.$peak_id, ]) %>%
  mutate(file_name=basename(fileNames(xdata_filled))[sample]) %>%
  select(-c(SNR, peak_cor, qscore))

split_groups <- split(feature_peaks, factor(feature_peaks$file_name))
files_newscores <- bplapply(seq_along(split_groups), 
                            function(x, split_xcms_filedata, ms_files, xdata,
                                     grabSingleFileData, xcmsQscoreCalculator,
                                     qscoreCalculator){
  library(data.table)
  file_peaks <- split_xcms_filedata[[x]]
  file_data <- grabSingleFileData(ms_files[x])
  file_data$rt <- xcms::adjustedRtime(xdata)[MSnbase::fromFile(xdata)==x][factor(file_data$rt)]
  file_data_table <- as.data.table(file_data)
  
  file_qscores <- lapply(seq_len(nrow(file_peaks)), FUN = xcmsQscoreCalculator, 
                         xcms_peakdf=file_peaks, file_data_table=file_data_table,
                         qscoreCalculator=qscoreCalculator)
  file_qscores_df <- do.call(rbind, file_qscores)
  return(cbind(split_xcms_filedata[[x]], file_qscores_df))
}, 
split_xcms_filedata=split_groups, 
ms_files=ms_files, 
xdata=xdata_filled,
grabSingleFileData = grabSingleFileData,
xcmsQscoreCalculator = xcmsQscoreCalculator,
qscoreCalculator = qscoreCalculator)
peaks_by_feature <- do.call(rbind, files_newscores) %>% 
  as.data.frame() %>% 
  filter(rt>60 & rt < 1100) %>%
  arrange(feature)
saveRDS(peaks_by_feature, file = "XCMS/temp_data/peaks_by_feature.rds")
print(Sys.time()-start_time)
# 2.5 minutes
beep(2)



### Remove isotope and adduct features ----
start_time <- Sys.time()
peaks_by_feature <- readRDS(file = "XCMS/temp_data/peaks_by_feature.rds")
xdata_filled <- readRDS(file = "XCMS/temp_data/current_xdata_filled.rds")

is_peak_iso <- bplapply(split(peaks_by_feature, peaks_by_feature$file_name), 
                        FUN = isIso, xdata=xdata_filled,
                        grabSingleFileData=grabSingleFileData,
                        checkPeakCor=checkPeakCor, pmppm=pmppm) %>%
  do.call(what = rbind) %>% as.data.frame()
is_peak_adduct <- bplapply(split(peaks_by_feature, peaks_by_feature$file_name), 
                           FUN = isAdduct, xdata=xdata_filled,
                           grabSingleFileData=grabSingleFileData,
                           checkPeakCor=checkPeakCor, pmppm=pmppm) %>%
  do.call(what = rbind) %>% as.data.frame()

addisod_peaks_by_feature <- peaks_by_feature %>%
  left_join(is_peak_iso) %>%
  left_join(is_peak_adduct)

addiso_features <- addisod_peaks_by_feature %>%
  group_by(feature) %>%
  summarise(prob_M1=mean(M1_match), prob_M2=mean(M2_match), prob_M3=mean(M3_match),
            prob_Na=mean(Na_match), prob_NH4=mean(NH4_match), 
            prob_H2O_H=mean(H2O_H_match), prob_2H=mean(`2H_match`)) %>%
  `[`(,-1) %>% `>`(0.8) %>% rowSums() %>% as.logical() %>% which() %>%
  `[`(unique(addisod_peaks_by_feature$feature), .)

cleaned_peaks_by_feature <- peaks_by_feature %>%
  filter(!feature%in%addiso_features)

removed_features <- addisod_peaks_by_feature %>%
  filter(feature%in%addiso_features) %>%
  group_by(feature) %>%
  summarise(mean_m_H_mz=mean(mz), mean_rt=mean(rt), mean_qscore=mean(qscore),
            prob_M1=mean(M1_match), prob_M2=mean(M2_match), prob_M3=mean(M3_match),
            prob_Na=mean(Na_match), prob_NH4=mean(NH4_match), 
            prob_H2O_H=mean(H2O_H_match), prob_2H=mean(`2H_match`))
removed_features[,-1] <- round(removed_features[,-1], digits = 4)
removed_features[removed_features<0.8] <- "-------"
print(removed_features)

saveRDS(cleaned_peaks_by_feature, file = "XCMS/temp_data/cleaned_peaks_by_feature.rds")
print(Sys.time()-start_time)
# 10 minutes
beep(2)



### Pause to plot all the "real" features! ----
library(plotly)
gp <- cleaned_peaks_by_feature %>%
  group_by(feature) %>%
  summarise(mean_mz=mean(mz), mean_rt=mean(rt),
            mean_rt_min=mean(rtmin), mean_rt_max=mean(rtmax),
            logmean_area=log10(mean(into))) %>%
  ggplot() + 
  geom_segment(aes(x=mean_rt_min, y = mean_mz, 
                   xend=mean_rt_max, yend=mean_mz,
                   color=logmean_area)) +
  scale_color_viridis_c()
ggplotly(gp)

#Run some preliminary ANOVA/Tukey (needs revision!)
# v <- cleaned_peaks_by_feature %>%
#   mutate(filenames = gsub("190715_", "", file_name)) %>%
#   left_join(metadata@data, by="filenames") %>%
#   filter(!depth%in%c("Pooled", "Std")) %>%
#   split(.$feature) %>%
#   lapply(FUN = function(x){
#     if(length(unique(x$depth))==1){
#       
#     }
#     TukeyHSD(aov(x$into~x$depth))$`x$depth`[,"p adj", drop=FALSE] %>%
#       cbind("depth.depth"=rownames(.)) %>%
#       merge(data.frame("p adj"=NA, `depth.depth`=c("Blank-DCM", "DCM=25m", "DCM-Blank")), 
#             by="depth.depth", all.y=TRUE)
#     return(tukeyout)
#   }) %>%
#   lapply(t)
  


### Calculate isotopes and adducts for remaining features ----
start_time <- Sys.time()
xdata_filled <- readRDS(file = "XCMS/temp_data/current_xdata_filled.rds")
cleaned_peaks_by_feature <- readRDS("XCMS/temp_data/cleaned_peaks_by_feature.rds")

all_file_peaks <- split(cleaned_peaks_by_feature, cleaned_peaks_by_feature$file_name)
peaks_isoadded <- bplapply(all_file_peaks, FUN = findIsoAdds, xdata=xdata_filled,
                           grabSingleFileData=grabSingleFileData, 
                           getPeakArea=getPeakArea, pmppm=pmppm, trapz=trapz) %>%
  do.call(what = rbind) %>% as.data.frame()

features_isoadded <- peaks_isoadded %>%
  group_by(feature) %>%
  summarise(mean_mz=mean(mz), mean_rt=mean(rt), area=mean(into), mean_q=mean(qscore),
            M1_area=mean(M_1_area)/area, prob_M1=mean(M_1_match), 
            M2_area=mean(M_2_area)/area, prob_M2=mean(M_2_match), 
            M3_area=mean(M_3_area)/area, prob_M3=mean(M_3_match),
            S34_area=mean(S_34_area)/area, prob_S34=mean(S_34_match),
            Na_area=mean(Na_area), prob_Na=mean(Na_match), 
            NH4_area=mean(NH4_area), prob_NH4=mean(NH4_match), 
            H2O_H_area=mean(H2O_H_area), prob_H2O_H=mean(H2O_H_match),
            X2H_area=mean(X2H_area), prob_2H=mean(X2H_match)) %>%
  #arrange(desc(rowSums(select(., starts_with("prob"))))) %>% 
  arrange(mean_mz) %>%
  as.data.frame()

write.csv(peaks_isoadded, file = "XCMS/temp_data/peaks_isoadded.csv", row.names = FALSE)
write.csv(features_isoadded, file = "XCMS/temp_data/features_isoadded.csv", row.names = FALSE)
print(Sys.time()-start_time)
# 10 minutes
beep(2)



### Assign molecular formulae ----
peak_df <- read.csv(file = "XCMS/temp_data/peaks_isoadded.csv", stringsAsFactors = FALSE)
feature_df <- read.csv(file = "XCMS/temp_data/features_isoadded.csv", stringsAsFactors = FALSE)
pubchem_formulas <- readRDS(file = "XCMS/unique_formulae.rds")
molecule_guesses <- pblapply(split(feature_df, feature_df$feature), function(x){
  mz <- x$mean_mz
  rdout <- Rdisop::decomposeMass(mz-1.007276, maxisotopes=3, ppm=5)
  if(is.null(rdout)){
    return(cbind(x, formula=NA, exactmass=NA))
  }
  rdout$isotopes <- NULL
  rdformat <- do.call(cbind, rdout) %>%
    as.data.frame() %>%
    select(formula, exactmass) %>%
    filter(sapply(formula, `%fin%`, pubchem_formulas))
  if(!nrow(rdformat)){
    return(cbind(x, formula=NA, exactmass=NA))
  }
  return(cbind(x, rdformat))
}) %>% do.call(what = rbind)
write.csv(molecule_guesses, file = "XCMS/temp_data/molecule_guesses.csv", row.names = FALSE)


# Plot # of 
molecule_guesses %>%
  group_by(feature, mean_mz, mean_rt) %>%
  summarise(num=sum(!is.na(formula))) %>%
  ggplot() + geom_jitter(aes(x=mean_mz, y=num), height=0.1) +
  theme_bw() + ylab("Number of possible formulae") + xlab("m/z")

molecule_guesses %>%
  group_by(feature, mean_mz, mean_rt) %>%
  summarise(num=sum(!is.na(formula))) %>%
  mutate(num=ifelse(num>1, 2, num)) %>%
  ggplot() + geom_jitter(aes(x=mean_rt, y=mean_mz, color=factor(num))) +
  theme_bw() + ylab("m/z") + xlab("Retention time") +
  scale_color_discrete(name="Number of\npossible\nformulae", labels=c(0, 1, "2+"))



### Add MS2 info ----
msms_files <- list.files("mzMLs/MSMS", pattern = "pos", full.names = TRUE)
raw_msmsdata <- pblapply(msms_files, grabSingleFileMS2)
nrgs <- as.numeric(gsub(".*neg|.*pos|\\.mzML", "", msms_files))
raw_msmsdata <- lapply(seq_along(raw_msmsdata), function(x){
  cbind(nrg=nrgs[x], raw_msmsdata[[x]])
})
raw_msmsdata <- as.data.table(do.call(rbind, raw_msmsdata))
feature_msms <- apply(feature_df, 1, FUN = function(feature_data){
  rtr <- as.numeric(feature_data["mean_rt"])+c(-10, 10)
  mzr <- pmppm(as.numeric(feature_data["mean_mz"]), ppm = 5)
  msms_data <- raw_msmsdata[rt%between%rtr & premz%between%mzr]
  if(!nrow(msms_data))return(NULL)
  return(cbind(feature=feature_data["feature"], msms_data))
})
v <- feature_msms %>%
  do.call(what = rbind) %>%
  left_join(feature_df, by="feature", all.y=TRUE) %>%
  left_join(molecule_guesses) %>%
  mutate(frag_found = apply(X = ., MARGIN = 1, FUN = fragCheck))

fragCheck <- function(row_data){
  formula <- row_data["formula"]
  frag_mz <- as.numeric(row_data["fragmz"])
  if(is.na(formula))return(NA)
  decomposeMass(mass = frag_mz, ppm = 5, maxElements = formula) %>%
    length() %>%
    as.logical()
}
  
frag_recovery <- pbsapply(seq_along(molecule_guesses$formula), function(i){
  if(is.na(molecule_guesses$formula[i]))return(0)
  if(!nrow(feature_msms[[i]]))return(0)
  frag_matches <- lapply(feature_msms[[i]]$fragmz, decomposeMass, ppm=5, 
                         maxElements=molecule_guesses$formula[i],
                         maxisotopes=1)
  frag_recovery <- mean(sapply(frag_matches, length)>0)
})



#Recover or validate betaine formula with MS2 fragments
formulas_to_check <- molecule_guesses[[5]]$formula
for(formula in formulas_to_check){
  frag_forms <- sapply(feature_msms[[5]]$`35`$fragmz-1.007276, 
                       Rdisop::decomposeMass, maxElements = formula, maxisotopes=1,
                       elements=gsub(pattern = "[[:digit:]]", replacement = "", formula),
                       ppm=5)
  print(frag_forms)
}
