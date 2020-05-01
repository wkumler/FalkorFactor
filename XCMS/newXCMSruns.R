# XCMS, but with CAMERA this time
# Just kidding, CAMERA doesn't do what I wanted it to at ALL
# Update: xMSannotator also sucks. Time to try SIRIUS!

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
# Make sure the dev version of XCMS is installed!




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
register(BPPARAM = SnowParam(tasks = length(ms_files), progressbar = TRUE))

metadata <- data.frame(
  fileid=basename(ms_files),
  sample_group=c("Blank", "Pooled", "Sample", "Std")[c(1, rep(2, 6), rep(3, 24), rep(4, 10))],
  depth=c("Blank", "Pooled", "DCM", "25m", "Std")[c(1, rep(2, 6), rep(c(rep(3, 3), rep(4, 3)), 4), rep(5, 10))],
  spindir=c("Blank", "Pooled", "Cyclone", "Anticyclone", "Std")[c(1, rep(2, 6), rep(3, 12), rep(4, 12), rep(5, 10))],
  time=c("Blank", "Pooled", "Morning", "Afternoon", "Std")[c(1, rep(2, 6), rep(c(rep(3, 6), rep(4, 6)), 2), rep(5, 10))]
) %>% new(Class = "NAnnotatedDataFrame")



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
saveRDS(xdata, file = "XCMS/temp_data/current_xdata.rds")
message(Sys.time()-start_time)
# 13 minutes

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
message(Sys.time()-start_time)
# 20 minutes



# Plot picked peaks, just for kicks! ----
xdata_cleanpeak <- readRDS("XCMS/temp_data/xdata_cleanpeak.rds")
given_file <- ms_files[2]
clean_peakdf <- chromPeaks(xdata_cleanpeak) %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  mutate(filename=fileNames(xdata_cleanpeak)[.[["sample"]]]) %>%
  arrange(mz) %>%
  filter(filename==given_file)
file_eic <- as.data.table(grabSingleFileData(given_file))
pdf(width = 17, height = 11, file = "XCMS/some_picked_peaks.pdf")
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
message(Sys.time()-start_time)
# 30 minutes

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
saveRDS(addiso_feature_defs, file = "XCMS/temp_data/isoadd_features.rds")
message(Sys.time()-start_time)
#40 minutes


# Calculate isotopes and adducts for remaining peaks ----
xdata_filled <- readRDS("XCMS/temp_data/current_xdata_filled.rds")
feature_defs <- featureDefinitions(xdata_filled)
addiso_feature_defs <- readRDS("XCMS/temp_data/isoadd_features.rds")

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
            O18_cor=median(O18_match), K_cor=median(K_match),
            S33_cor=median(S33_match),
            Na_cor=median(Na_match), NH4_cor=median(NH4_match), 
            H2O_H_cor=median(H2O_H_match), X2H_cor=median(X2H_match)) %>% 
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
  select(1, grep("R2", names(peakareamatch))) %>%
  pivot_longer(cols = -c("feature"), names_to = "addiso", values_to = "R2") %>%
  mutate(addiso=gsub("_R2", "", addiso))
peak_slopes <- peak_slope_R2 %>%
  select(1, grep("slope", names(peakareamatch))) %>%
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
  mutate(rel_int=ifelse(cor>0.8&R2>0.9, slope*avgarea, 0)) %>%
  select(-c("cor", "R2", "slope")) %>%
  pivot_wider(names_from = addiso, values_from = rel_int)

saveRDS(final_peaks, file = "XCMS/final_peaks.rds")
saveRDS(final_features, file = "XCMS/final_features.rds")


# Analysis! ----
# Single compound annotation sanity checks
final_peaks <- readRDS(file = "XCMS/final_peaks.rds")

# Formula: C5H12NO2 (betaine+H)
feature_num <- "FT054"
options(scipen = 5)

ft_isodata <- final_peaks %>% filter(feature==feature_num) %>%
  select(M_area, C13_area, N15_area, O18_area, X2C13_area, S34_area, S33_area) %>%
  pivot_longer(cols = starts_with(c("C", "X", "N", "O", "S")))

ggplot(ft_isodata, aes(x=M_area, y=value)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_wrap(~name, scales = "free_y") +
  ggtitle(feature_num)

lmoutput <- split(ft_isodata, ft_isodata$name) %>%
  lapply(lm, formula=value~M_area) %>%
  lapply(summary) %>%
  lapply(`[[`, "coefficients")

(0:10)[which.min(abs(sapply(0:10, dbinom, x=1, prob=0.0107)-lmoutput[["C13_area"]]["M_area", "Estimate"]))]
(0:10)[which.min(abs(sapply(0:10, dbinom, x=2, prob=0.0107)-lmoutput[["X2C13_area"]]["M_area", "Estimate"]))]
(0:10)[which.min(abs(sapply(0:10, dbinom, x=1, prob=0.00368)-lmoutput[["N15_area"]]["M_area", "Estimate"]))]
(0:10)[which.min(abs(sapply(0:10, dbinom, x=1, prob=0.00205)-lmoutput[["O18_area"]]["M_area", "Estimate"]))]
(0:10)[which.min(abs(sapply(0:10, dbinom, x=1, prob=0.0421)-lmoutput[["S34_area"]]["M_area", "Estimate"]))]
(0:10)[which.min(abs(sapply(0:10, dbinom, x=1, prob=0.0075)-lmoutput[["S33_area"]]["M_area", "Estimate"]))]
final_peaks %>% filter(feature==feature_num) %>%
  summarize(mzmed=median(mz), rtmed=median(rt)) %>%
  mutate(C13=mzmed+1.003355, X2C13=mzmed+1.003355*2,
         N15=mzmed+0.997035, O18=mzmed+2.004244, 
         S33=mzmed+0.999387, S34=mzmed+1.995796)


seq(0.005, 0.015, 0.0001)[which.min(abs(sapply(seq(0.005, 0.015, 0.0001), dbinom, x=1, size=5)-lmoutput[["C13_area"]]["M_area", "Estimate"]))]


# Multi-peak stuff
final_peaks <- readRDS(file = "XCMS/final_peaks.rds")
final_features <- final_peaks %>% 
  group_by(feature) %>%
  summarize(mzmed=median(mz), rtmed=median(rt), avgarea=mean(M_area)) %>%
  as.data.frame(stringsAsFactors=FALSE)
final_diffreport <- split(final_peaks, final_peaks$feature) %>%
  lapply(function(x){
    DCM_areas <- x$M_area[grep(pattern = "DCM", x$file_name)]
    m25_areas <- x$M_area[grep(pattern = "25m", x$file_name)]
    if(length(DCM_areas)<3|length(m25_areas)<3)return(c(0.05, 1))
    c(pval=t.test(DCM_areas, m25_areas)$p.value, 
      diff=mean(DCM_areas)/mean(m25_areas))
  }) %>% do.call(what = rbind) %>% as.data.frame() %>%
  mutate(feature=rownames(.)) %>% left_join(final_features, ., by="feature") %>%
  arrange(mzmed)
DCM_enriched <- final_diffreport %>%
  filter(diff>1) %>%
  filter(pval<0.01) %>%
  arrange(pval)
surface_enriched <- final_diffreport %>%
  filter(diff<1) %>%
  filter(pval<0.01) %>%
  arrange(pval)



# SIRIUS ----
final_peaks <- readRDS(file = "XCMS/final_peaks.rds")
final_features <- final_peaks %>% 
  group_by(feature) %>%
  summarize(mzmed=median(mz), rtmed=median(rt), avgarea=mean(M_area),
            areamed=median(M_area), 
            C13_areamed=ifelse(median(C13_match>0.9), median(C13_area), NA), 
            X2C13_areamed=ifelse(median(X2C13_match>0.9), median(X2C13_area), NA), 
            S34_areamed=ifelse(median(S34_match>0.9), median(S34_area), NA),
            S33_areamed=ifelse(median(S33_match>0.9), median(S33_area), NA),
            N15_areamed=ifelse(median(N15_match>0.9), median(N15_area), NA),
            O18_areamed=ifelse(median(O18_match>0.9), median(O18_area), NA),
            Na_areamed=ifelse(median(Na_match>0.9), median(Na_area), NA),
            NH4_areamed=ifelse(median(NH4_match>0.9), median(NH4_area), NA),
            K_areamed=ifelse(median(K_match>0.9), median(K_area), NA),
            H2O_H_areamed=ifelse(median(H2O_H_match>0.9), median(H2O_H_area), NA),
            X2H_areamed=ifelse(median(X2H_match>0.9), median(X2H_area), NA)) %>%
  as.data.frame(stringsAsFactors=FALSE)

MSMS_files <- "mzMLs/MSMS/" %>%
  list.files(pattern = ".mzML", full.names = TRUE) %>%
  normalizePath() %>%
  `[`(grepl("pos", x = .))
raw_MSMS_data <- lapply(MSMS_files, grabSingleFileMS2) %>%
  mapply(FUN = cbind, as.numeric(gsub(".*DDApos|.mzML", "", MSMS_files)), 
         SIMPLIFY = FALSE) %>%
  lapply(`names<-`, c("rt", "premz", "fragmz", "int", "voltage")) %>%
  do.call(what = rbind) %>% as.data.table()


# MGF method
output_dir <- "XCMS/sirius_temp"
mgf_maker <- function(feature_msdata, ms1, ms2, output_file){
  if(!nrow(ms2)){
    outtext <- c("BEGIN IONS",
                 paste0("PEPMASS=", feature_msdata$mzmed),
                 "MSLEVEL=1",
                 "CHARGE=1+",
                 apply(ms1, 1, paste, collapse=" "),
                 "END IONS",
                 "")
  } else {
    outtext <- c("BEGIN IONS",
                 paste0("PEPMASS=", feature_msdata$mzmed),
                 "MSLEVEL=1",
                 "CHARGE=1+",
                 apply(ms1, 1, paste, collapse=" "),
                 "END IONS",
                 "",
                 "BEGIN IONS",
                 paste0("PEPMASS=", feature_msdata$mzmed),
                 "MSLEVEL=2",
                 "CHARGE=1+",
                 apply(ms2[ms2$voltage==20, c("fragmz", "int")], 
                       1, paste, collapse=" "),
                 "END IONS",
                 "",
                 "BEGIN IONS",
                 paste0("PEPMASS=", feature_msdata$mzmed),
                 "MSLEVEL=2",
                 "CHARGE=1+",
                 apply(ms2[ms2$voltage==35, c("fragmz", "int")], 
                       1, paste, collapse=" "),
                 "END IONS",
                 "BEGIN IONS",
                 paste0("PEPMASS=", feature_msdata$mzmed),
                 "MSLEVEL=2",
                 "CHARGE=1+",
                 apply(ms2[ms2$voltage==50, c("fragmz", "int")], 
                       1, paste, collapse=" "),
                 "END IONS")
  }
  writeLines(outtext, con = output_file)
}

for(feature_num in head(final_features$feature, 40)){
  output_file <- paste0(output_dir, "\\", feature_num, ".mgf")
  feature_msdata <- final_features[final_features$feature==feature_num, ]
  ms1 <- rbind(c(feature_msdata$mzmed, feature_msdata$areamed),
               c(feature_msdata$mzmed+1.003355, feature_msdata$M1_areamed),
               c(feature_msdata$mzmed+1.003355*2, feature_msdata$M2_areamed))
  ms1 <- ms1[!is.na(ms1[,2]), , drop=FALSE]
  ms2 <- raw_MSMS_data[premz%between%pmppm(feature_msdata$mzmed)&
                         rt%between%(feature_msdata$rtmed+c(-10, 10))]
  mgf_maker(feature_msdata = feature_msdata, ms1 = ms1, 
            ms2 = ms2, output_file = output_file)
}


sirius_cmd <- paste0('"C://Program Files//sirius-win64-4.0.1//',
                     'sirius-console-64.exe" ',
                     ' -p orbitrap',
                     ' --database bio',
                     #' --fingerid',
                     ' -i [M+H]+ ', '"', normalizePath(output_dir), '"')

sirius_output <- system(sirius_cmd, intern = TRUE)
sirius_clean <- sirius_output %>%
  grep(pattern = "^Sirius|^[[:digit:]]", value = TRUE) %>%
  split(cumsum(grepl(pattern = "Sirius", x = .))) %>%
  lapply(function(feature_data){
    dput(feature_data)
    feature_num <- substr(x = feature_data[1], start = 22, 26)
    clean_data <- lapply(strsplit(feature_data[-1], split = "\t"), function(x){
      x[1] <- sub(pattern = "^[[:digit:]]\\.\\)\ ", "", x = x[1])
      x <- gsub(pattern = ".*\\:\\ \\+*", "", x)
      x <- gsub(pattern = "\\ \\%", "", x)
    })
    if(!length(clean_data)){
      output <- matrix(c(feature_num, rep(NA, 8)), nrow = 1)
    } else {
      output <- do.call(rbind, clean_data) %>%
        cbind(feature_num, .)
    }
    colnames(output) <- c("feature", "formula", "adduct", "score", "tree", 
                          "iso", "peaks", "expl_int", "iso_peaks")
    return(output)
  }) %>% do.call(what = rbind)
sirius_clean

sirius_clean %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  left_join(final_features) %>%
  arrange(mzmed) %>% select(!contains("_"))

final_features %>% filter(feature=="FT002")
