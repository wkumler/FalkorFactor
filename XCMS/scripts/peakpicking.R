# Script to run XCMS peakpicking, retention time correction, and grouping
# followed by de-adducting/de-isotoping
# Custom de-adducting and de-isotoping because CAMERA and xMSannotator suck
# Could be improved by including a triple-check via clustering function
# Expects .mzML files in mzMLs subdirectory of project dir
# Returns a de-adducted/de-isotoped peak list
# Requires dev version of XCMS for improved peakpicking given settings



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
  file_peaks <- file_peaks[order(file_peaks$rtmax),]
  file_path <- paste("mzMLs", unique(file_peaks$file_name), sep = "/")
  file_data <- grabSingleFileData(file_path)
  file_data$rt <- xcms::adjustedRtime(xdata)[
    MSnbase::fromFile(xdata)==unique(file_peaks$sample)][factor(file_data$rt)]
  library(data.table)
  file_data_table <- as.data.table(file_data)
  
  peak_splits <- split(file_peaks, ceiling(seq_len(nrow(file_peaks))/10))
  
  iso_matches_all <- lapply(peak_splits, function(i){
    eic_many <- file_data_table[rt>min(i$rtmin)&rt<max(i$rtmax)]
    individual_peaks <- split(i, seq_len(nrow(i)))
    iso_matches <- lapply(individual_peaks, function(peak_row_data){
      init_eic <- eic_many[mz%between%pmppm(peak_row_data["mz"], ppm = 5) & 
                             rt%between%c(peak_row_data["rtmin"], peak_row_data["rtmax"])]
      init_area <- trapz(init_eic$rt, init_eic$int)
      
      
      isos_to_check <- c(C13=-1.003355, X2C13=-2*1.003355, S34=-1.995796, S33=-0.999387,
                         N15=-0.997035, O18=-2.004244) + peak_row_data[["mz"]]
      adducts_to_check <- c(Na=-22.98922+1.007276, NH4=-18.0338+1.007276,
                            H2O_H=+18.0106, K=-38.963708+1.007276) + peak_row_data[["mz"]]
      more_adducts <- c(X2H=peak_row_data[["mz"]]*2-2*1.007276+1.007276)
      
      masses_to_check <- c(isos_to_check, adducts_to_check, more_adducts)
      
      output <- lapply(masses_to_check, checkPeakCor, rtmin=peak_row_data["rtmin"], 
                       rtmax=peak_row_data["rtmax"], init_eic = init_eic, 
                       file_data_table = eic_many, pmppm = pmppm, trapz = trapz)
      linear_output <- do.call(c, output)
      names(linear_output) <- paste0(rep(names(masses_to_check), each=2), 
                                     c("_match", "_area"))
      linear_output <- c(M_area=init_area, linear_output)
      return(linear_output)
    })
    iso_matches <- do.call(rbind, iso_matches)
    return(iso_matches)
  })
  iso_matches_all <- do.call(rbind, iso_matches_all)
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
findIsoAdduct <- function(file_peaks, xdata, grabSingleFileData,
                          checkPeakCor, pmppm, trapz){
  file_peaks <- file_peaks[order(file_peaks$rtmax),]
  file_path <- paste("mzMLs", unique(file_peaks$file_name), sep = "/")
  file_data <- grabSingleFileData(file_path)
  file_data$rt <- xcms::adjustedRtime(xdata)[
    MSnbase::fromFile(xdata)==unique(file_peaks$sample)][factor(file_data$rt)]
  library(data.table)
  file_data_table <- as.data.table(file_data)
  
  peak_splits <- split(file_peaks, ceiling(seq_len(nrow(file_peaks))/10))
  
  iso_matches_all <- lapply(peak_splits, function(i){
    eic_many <- file_data_table[rt>min(i$rtmin)&rt<max(i$rtmax)]
    individual_peaks <- split(i, seq_len(nrow(i)))
    iso_matches <- lapply(individual_peaks, function(peak_row_data){
      init_eic <- eic_many[mz%between%pmppm(peak_row_data["mz"], ppm = 5) & 
                             rt%between%c(peak_row_data["rtmin"], peak_row_data["rtmax"])]
      init_area <- trapz(init_eic$rt, init_eic$int)
      
      
      isos_to_check <- c(C13=1.003355, X2C13=2*1.003355, S34=1.995796, S33=0.999387,
                         N15=0.997035, O18=2.004244) + peak_row_data[["mz"]]
      adducts_to_check <- c(Na=22.98922-1.007276, NH4=18.0338-1.007276,
                            H2O_H=-18.0106, K=38.963708-1.007276) + peak_row_data[["mz"]]
      more_adducts <- c(X2H=peak_row_data[["mz"]]-1.007276+2*1.007276)/2
      
      masses_to_check <- c(isos_to_check, adducts_to_check, more_adducts)
      
      output <- lapply(masses_to_check, checkPeakCor, rtmin=peak_row_data["rtmin"], 
                       rtmax=peak_row_data["rtmax"], init_eic = init_eic, 
                       file_data_table = eic_many, pmppm = pmppm, trapz = trapz)
      linear_output <- do.call(c, output)
      names(linear_output) <- paste0(rep(names(masses_to_check), each=2), 
                                     c("_match", "_area"))
      linear_output <- c(M_area=init_area, linear_output)
      return(linear_output)
    })
    iso_matches <- do.call(rbind, iso_matches)
    return(iso_matches)
  })
  iso_matches_all <- do.call(rbind, iso_matches_all)
  return(cbind(file_peaks, iso_matches_all))
}



# Setup things ----
start_time <- Sys.time()
library(tidyverse)
library(data.table)
library(pbapply)
library(xcms)

ms_files <- "mzMLs" %>%
  list.files(pattern = ".mzML", full.names = TRUE) %>%
  normalizePath() %>%
  `[`(!grepl("Fullneg|Fullpos|QC-KM1906", x = .))

metadata <- data.frame(
  fileid=basename(ms_files),
  sample_group=c("Blank", "Pooled", "Sample", "Std")[c(1, rep(2, 6), rep(3, 24), rep(4, 10))],
  depth=c("Blank", "Pooled", "DCM", "25m", "Std")[c(1, rep(2, 6), rep(c(rep(3, 3), rep(4, 3)), 4), rep(5, 10))],
  spindir=c("Blank", "Pooled", "Cyclone", "Anticyclone", "Std")[c(1, rep(2, 6), rep(3, 12), rep(4, 12), rep(5, 10))],
  time=c("Blank", "Pooled", "Morning", "Afternoon", "Std")[c(1, rep(2, 6), rep(c(rep(3, 6), rep(4, 6)), 2), rep(5, 10))]
) %>% new(Class = "NAnnotatedDataFrame")
write.csv(x = metadata@data, file = "XCMS/data_pretty/metadata.csv", row.names = FALSE)

sirius_project_dir <- "XCMS/data_intermediate/sirius_project"
register(BPPARAM = SnowParam(tasks = length(ms_files), progressbar = TRUE))


# Peakpicking ----
raw_data <- readMSData(files = ms_files, pdata = metadata, msLevel. = 1, 
                       verbose = TRUE, centroided. = TRUE, mode = "onDisk")
cwp <- CentWaveParam(ppm = 2.5, peakwidth = c(15, 15), 
                     snthresh = 1, prefilter = c(0, 10000), 
                     integrate = 2, mzCenterFun = "wMean", 
                     mzdiff = 0.001, fitgauss = FALSE, 
                     noise = 5000, firstBaselineCheck = FALSE, 
                     extendLengthMSW = TRUE)
xdata <- suppressMessages(findChromPeaks(raw_data, param = cwp))
saveRDS(xdata, file = "XCMS/data_intermediate/current_xdata.rds")
message(Sys.time()-start_time)
# 13 minutes

xdata <- readRDS(file = "XCMS/data_intermediate/current_xdata.rds")
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
write.csv(x = peakdf_qscored, row.names = FALSE,
          file = "XCMS/data_intermediate/all_peaks_w_qscores.csv")
xdata_cleanpeak <- `chromPeaks<-`(xdata, peakdf_qscored)
saveRDS(xdata_cleanpeak, file = "XCMS/data_intermediate/xdata_cleanpeak.rds")
message(Sys.time()-start_time)
# 20 minutes



# Other XCMS things (rtcor, group) ----
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

saveRDS(xdata_filled, file = "XCMS/data_intermediate/current_xdata_filled.rds")
message(Sys.time()-start_time)
# 30 minutes

show(featureDefinitions(xdata_filled))



# Find isotopes and adducts ----
xdata_filled <- readRDS("XCMS/data_intermediate/current_xdata_filled.rds")
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
#6.5 minutes

peakshapematch <- is_peak_iso %>%
  group_by(feature) %>%
  summarise(C13_prob=median(C13_match), X2C13_prob=median(X2C13_match), 
            S34_prob=median(S34_match), S33_prob=median(S33_match),
            N15_prob=median(N15_match), O18_prob=median(O18_match),
            Na_prob=median(Na_match), NH4_prob=median(NH4_match), 
            H2O_H_prob=median(H2O_H_match), K_prob=median(K_match),
            X2H_prob=median(X2H_match)) %>%
  as.data.frame()

peakareamatch <- lapply(unique(is_peak_iso$feature), function(i){
  feature_areas <- is_peak_iso[is_peak_iso$feature==i,]
  area_cols <- grep(pattern = "area$", names(feature_areas), value = TRUE)[-1]
  sapply(area_cols, function(x){
    suppressWarnings(cor(feature_areas$M_area, feature_areas[[x]]))
  })
}) %>% do.call(what=rbind) %>% 
  `[<-`(is.na(.), 0) %>% 
  as.data.frame(stringsAsFactors=FALSE) %>%
  mutate(feature=unique(is_peak_iso$feature)) %>%
  select(feature, everything()) %>%
  arrange(feature)

likely_addisos <- peakareamatch$feature[
  which(rowSums(peakareamatch[,names(peakareamatch)!="feature"]>0.95&
                  peakshapematch[,names(peakshapematch)!="feature"]>0.9)>=1)
  ]

addiso_feature_defs <- feature_defs %>%
  `[`(peakareamatch$feature%in%likely_addisos, c("mzmed", "rtmed")) %>%
  as.data.frame() %>%
  round(digits = 5)

v <- peakareamatch[peakareamatch$feature%in%likely_addisos,]
v <- cbind(v, peakshapematch[peakshapematch$feature%in%likely_addisos,-1])
cmpds <- c("C13", "X2C13", "N15", "O18", "S33", "S34", "Na", "NH4", "K", "H2O_H", "X2H")
v <- v[,paste0(rep(cmpds, each=2), c("_prob", "_area"))]
v <- round(v, digits = 5)
v[v<0.9] <- "-------"
v <- cbind(addiso_feature_defs, v[,-1])

v[which(v["S33_area"]!="-------"&v["S33_prob"]!="-------"&v["rtmed"]<780),]
write.csv(addiso_feature_defs, row.names = FALSE,
          file = "XCMS/data_intermediate/isotope_adduct_features.csv")
message(Sys.time()-start_time)
#40 minutes


# Calculate isotopes and adducts for remaining peaks ----
xdata_filled <- readRDS("XCMS/data_intermediate/current_xdata_filled.rds")
feature_defs <- featureDefinitions(xdata_filled)
addiso_feature_defs <- read.csv("XCMS/data_intermediate/isotope_adduct_features.csv")

# Find the chromPeaks associated with each featureDefinition
# Remove the features identified in addiso_feature_defs as likely adducts/isotopes
clean_feature_peaks <- lapply(seq_len(nrow(feature_defs)), function(i){
  cbind(feature=sprintf("FT%03d", i), 
        peak_id=unlist(feature_defs$peakidx[i]))
}) %>% 
  do.call(what=rbind) %>% 
  as.data.frame(stringsAsFactors=FALSE) %>% 
  mutate(peak_id=as.numeric(peak_id)) %>%
  cbind(chromPeaks(xdata_filled)[.$peak_id, ]) %>%
  mutate(file_name=basename(fileNames(xdata_filled))[sample]) %>%
  arrange(feature, sample) %>%
  filter(!feature%in%rownames(addiso_feature_defs))

# For each peak, look for data at +/- each adduct/isotope m/z 
# Also calculate cor while the raw data is being accessed anyway
split_list <- split(clean_feature_peaks, clean_feature_peaks$file_name)
final_peaks <- bplapply(split_list, FUN = findIsoAdduct, xdata=xdata_filled,
                        grabSingleFileData=grabSingleFileData,
                        checkPeakCor=checkPeakCor, 
                        pmppm=pmppm, trapz=trapz) %>%
  do.call(what = rbind) %>% as.data.frame() %>% 
  `rownames<-`(NULL) %>% arrange(feature)

# Calculate median cor for each FEATURE from the various peak cors
peak_cors <- final_peaks %>%
  group_by(feature) %>%
  summarise(C13_cor=median(C13_match), X2C13_cor=median(X2C13_match), 
            S34_cor=median(S34_match), N15_cor=median(N15_match), 
            O18_cor=median(O18_match), S33_cor=median(S33_match),
            K_cor=median(K_match), Na_cor=median(Na_match), 
            NH4_cor=median(NH4_match), H2O_H_cor=median(H2O_H_match), 
            X2H_cor=median(X2H_match)) %>% 
  pivot_longer(cols = -c("feature"), names_to = "addiso", values_to = "cor") %>%
  mutate(addiso=gsub("_cor", "", addiso))
# For each feature, plot adduct/iso areas against OG peak areas
# Run lm() to get best fit line slope and R-squared
peak_slope_R2 <- lapply(unique(final_peaks$feature), function(i){
  feature_areas <- final_peaks[final_peaks$feature==i,]
  area_cols <- grep(pattern = "area$", names(feature_areas), value = TRUE)[-1]
  area_outputs <- lapply(area_cols, function(x){
    lmoutput <- lm(feature_areas[[x]]~feature_areas$M_area)
    useful_info <- c(r2=summary(lmoutput)$r.squared, 
                     slope=lmoutput$coefficients["feature_areas$M_area"])
    return(useful_info)
  })
  lapply(area_outputs, c) %>%
    unlist() %>%
    `names<-`(paste0(rep(gsub("area", "", area_cols), each=2), c("R2", "slope")))
}) %>% 
  do.call(what=rbind) %>% `[<-`(is.na(.), 0) %>% 
  as.data.frame(stringsAsFactors=FALSE) %>%
  mutate(feature=unique(final_peaks$feature)) %>%
  select(feature, everything()) %>%
  arrange(feature)
# Separate out R-squareds and slopes (easier to do here than after merging)
peak_R2s <- peak_slope_R2 %>%
  select(1, grep("R2", names(peak_slope_R2))) %>%
  pivot_longer(cols = -c("feature"), names_to = "addiso", values_to = "R2") %>%
  mutate(addiso=gsub("_R2", "", addiso))
peak_slopes <- peak_slope_R2 %>%
  select(1, grep("slope", names(peak_slope_R2))) %>%
  pivot_longer(cols = -c("feature"), names_to = "addiso", values_to = "slope") %>%
  mutate(addiso=gsub("_slope", "", addiso))


# Establish thresholds for "yes, this is probably an adduct"
# If above threshold, return peak area as relative intensity
# If below, return nothing
# Essentially produces a cleaned up MS1 spectrum with only adducts/isotopes
final_features <- final_peaks %>% 
  group_by(feature) %>%
  summarize(mzmed=median(mz), rtmed=median(rt), avgarea=mean(M_area)) %>%
  left_join(peak_cors, by="feature") %>%
  left_join(peak_R2s, by=c("feature", "addiso")) %>%
  left_join(peak_slopes, by=c("feature", "addiso")) %>%
  mutate(rel_int=ifelse(cor>0.8&R2>0.9, round(slope*avgarea), 0)) %>%
  select(-c("cor", "R2", "slope")) %>%
  pivot_wider(names_from = addiso, values_from = rel_int)



# Write out peak and feature lists ----
write.csv(x = final_peaks, file = "XCMS/data_pretty/final_peaks.csv")
write.csv(x = final_features, file = "XCMS/data_pretty/final_features.csv")


