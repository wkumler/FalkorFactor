
# Setup things ----
library(xcms)
library(tidyverse)
library(data.table)
library(beepr)
library(httr)
library(pbapply)
start_time <- Sys.time()

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
) %>% new(Class = "NAnnotatedDataFrame")



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
qscoreCalculator <- function(eic){
  #Check for bogus EICs
  if(nrow(eic)<5){
    return(data.frame(SNR=0, peak_cor=0, qscore=0))
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
  output <- data.frame(SNR, peak_cor, qscore=SNR*peak_cor^2*log10(max(eic$int)))
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
searchPubchem <- function(mass, ppm, allowed_atoms=c("C", "H", "N", "O", "P", "S")){
  baseurl <- "https://pubchem.cheminfo.org/mfs/em?em=%f&precision=%.2f"
  url <- sprintf(baseurl, mass, ppm)
  response <- GET(url)
  if(!response$status_code==200){
    stop(paste0("Trouble accessing PubChem, error code:", response$status_code))
  }
  raw_content <- content(response)
  removed_compounds <- 0
  clean_content <- lapply(raw_content$result, function(x, allowed_atoms){
    atoms <- as.data.frame(x$atom)
    if(any(!names(atoms)%in%allowed_atoms)){
      removed_compounds <<- removed_compounds+1
      return(NULL)
    }
    x$atom <- NULL
    x$unsaturation <- NULL
    as.data.frame(x, stringsAsFactors=FALSE)
  }, allowed_atoms=allowed_atoms)
  content_df <- do.call(rbind, clean_content)
  print(paste("Found", nrow(content_df), "reasonable formula"))
  if(removed_compounds)print(paste("Removed", removed_compounds, "with weird atoms"))
  rm(removed_compounds)
  return(content_df)
}
findAdducts <- function(file_peaks, xdata){
  file_path <- paste("mzMLs", unique(file_peaks$file_name), sep = "/")
  file_data <- grabSingleFileData(file_path)
  file_data$rt <- xcms::adjustedRtime(xdata)[
    MSnbase::fromFile(xdata)==unique(file_peaks$sample)][
      factor(file_data$rt)]
  file_data_table <- as.data.table(file_data)
  
  outlist <- list()
  for(i in seq_len(nrow(file_peaks))){
    peak_row_data <- file_peaks[i, ]
    mzrange <- c(peak_row_data$mzmin, peak_row_data$mzmax)
    rtrange <- c(peak_row_data$rtmin, peak_row_data$rtmax)
    init_eic <- file_data_table[mz%between%mzrange & rt%between%rtrange]
    is_M1_iso <- checkForPeak(mass = peak_row_data$mz-1.003355, 
                              rtmin = peak_row_data$rtmin, 
                              rtmax = peak_row_data$rtmax, 
                              init_eic = init_eic, 
                              file_data=file_data_table)
    is_sodium <- checkForPeak(mass = peak_row_data$mz-21.98249, 
                              rtmin = peak_row_data$rtmin, 
                              rtmax = peak_row_data$rtmax, 
                              init_eic = init_eic, 
                              file_data=file_data_table)
    outlist[[i]] <- data.frame(is_M1_iso, is_sodium)
  }
  return(do.call(rbind, outlist))
}
checkForPeak <- function(mass, rtmin, rtmax, init_eic, file_data,
                         quality_cutoff = 1, similarity_cutoff=0.8){
  given_eic <- file_data[mz%between%pmppm(mass, ppm = 5) & rt%between%c(rtmin, rtmax)]
  if(nrow(given_eic)<3){
    return(FALSE)
  }
  peak_qscore <- qscoreCalculator(given_eic)$qscore
  
  merged_eic <- merge(init_eic, given_eic, by="rt")
  if(nrow(merged_eic)<3){
    return(FALSE)
  }
  peak_match <- cor(merged_eic$int.x, merged_eic$int.y)
  
  if(peak_qscore>quality_cutoff & peak_match>similarity_cutoff){
    return(TRUE)
  } else {
    return(FALSE)
  }
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
# Load ms_files here for debug purposes
ms_files <- "mzMLs" %>%
  list.files(pattern = ".mzML", full.names = TRUE) %>%
  normalizePath() %>%
  `[`(!grepl("Fullneg|Fullpos|QC-KM1906", x = .))



# Load MS data ----
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



# Perform peakpicking ----
start_time <- Sys.time()
raw_data <- readRDS(file = "XCMS/temp_data/current_raw_data.rds")
register(BPPARAM = SnowParam(tasks = length(ms_files), progressbar = TRUE))
cwp <- CentWaveParam(ppm = 5, peakwidth = c(20, 80), 
                     snthresh = 0, prefilter = c(0, 0), 
                     integrate = 1, mzCenterFun = "wMean", 
                     mzdiff = 0.0001, fitgauss = FALSE, 
                     noise = 0, firstBaselineCheck = FALSE)
xdata <- suppressMessages(findChromPeaks(raw_data, param = cwp))
print(xdata)
saveRDS(xdata, file = "XCMS/temp_data/current_xdata.rds")
print(Sys.time()-start_time)
# 30 minutes
beep(2)



# Re-assign quality scores to confirm good peaks ----
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



# Decide on quality threshold and reassign peaklist ----
start_time <- Sys.time()
xdata <- readRDS(file = "XCMS/temp_data/current_xdata.rds")
peakdf_qscored <- as.matrix(read.csv(file = "XCMS/temp_data/peakdf_qscored.csv"))
threshold <- 100
cleandf_qscored <- peakdf_qscored[peakdf_qscored[, "qscore"]>threshold,]
xdata_cleanpeak <- `chromPeaks<-`(xdata, cleandf_qscored)
# sample_df <- peakdf_qscored[sample(1:nrow(peakdf_qscored), size = 10000),]
# plot(log10(sample_df[,"sn"]), log10(sample_df[,"qscore"]))



# Adjust retention time and compare ----
start_time <- Sys.time()
obp <- ObiwarpParam(binSize = 0.01, centerSample = 4, response = 1, distFun = "cor_opt")
xdata_rt <- suppressMessages(adjustRtime(xdata_cleanpeak, param = obp))
plotAdjustedRtime(xdata_rt)
saveRDS(xdata_rt, file = "XCMS/temp_data/current_xdata_rt.rds")
print(Sys.time()-start_time)
# 16 minutes
beep(2)



# Correspondence ----
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



# Fill peaks ----
start_time <- Sys.time()
xdata_cor <- readRDS(file = "XCMS/temp_data/current_xdata_cor.rds")
xdata_filled <- suppressMessages(fillChromPeaks(xdata_cor, param = FillChromPeaksParam()))
saveRDS(object = xdata_filled, file = "XCMS/temp_data/current_xdata_filled.rds")
# 1 minute
print(Sys.time()-start_time)


# Add new quality scores ----
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
feature_peaks_rescored <- do.call(rbind, files_newscores) %>% 
  as.data.frame() %>% arrange(feature)
saveRDS(feature_peaks_rescored, file = "XCMS/temp_data/feature_peaks_rescored.rds")
print(Sys.time()-start_time)
# 1 minute
beep(2)

# Identify and remove adducts and isotope features ----
start_time <- Sys.time()
feature_peaks_rescored <- readRDS("XCMS/temp_data/feature_peaks_rescored.rds")
xdata_filled <- readRDS(file = "XCMS/temp_data/current_xdata_filled.rds")
all_file_peaks <- split(feature_peaks_rescored, feature_peaks_rescored$file_name)
feature_adduct_iso <- pblapply(all_file_peaks, FUN = findAdducts, xdata=xdata_filled)
# Takes ~3 minutes
feature_adduct_iso <- do.call(rbind, feature_adduct_iso)
feature_peaks_added <- feature_peaks_rescored %>% 
  left_join(feature_adduct_iso, by="peak_id") %>%
  split(.$feature) %>%
  lapply(FUN = merge, 
         y=data.frame(file_name=basename(ms_files)), 
         by=c("file_name"), all.y=TRUE) %>%
  do.call(what = rbind) %>%
  mutate(feature=rep(unique(feature_peaks_rescored$feature), each=length(ms_files)))
saveRDS(feature_peaks_added, file = "XCMS/temp_data/feature_peaks_added.rds")
print(Sys.time()-start_time)
# 2.5 minutes
beep(2)

feature_peaks_added <- readRDS("XCMS/temp_data/feature_peaks_added.rds")
zscore <- function(blank, samples){
  sig <- sd(samples, na.rm = TRUE)/sqrt(sum(!is.na(samples)))
  zscore <- (mean(samples, na.rm=TRUE)-blank)/sig
  return(zscore)
}
feature_summaries <- feature_peaks_added %>%
  `[<-`(is.na(.), value=0) %>%
  group_by(feature) %>%
  summarise(mz=weighted.mean(mz, qscore), rt=weighted.mean(rt, qscore),
            M1_iso=weighted.mean(M_1_area/into, M_1_q*M_1_simil^2, na.rm = TRUE),
            M2_iso=weighted.mean(M_2_area/into, M_2_q*M_2_simil^2, na.rm = TRUE)) %>%
  as.data.frame()

plotPeak <- function(peakdf_row){
  v <- ms_files[feature_peaks_added[peakdf_row,"sample"]] %>%
    grabSingleFileData() %>%
    mutate(rt = xcms::adjustedRtime(xdata_filled)[
      MSnbase::fromFile(xdata_filled)==feature_peaks_added[peakdf_row,"sample"]
      ][factor(rt)]) %>%
    filter(rt>feature_peaks_added[peakdf_row,"rtmin"]&
             rt<feature_peaks_added[peakdf_row,"rtmax"]&
             mz>feature_peaks_added[peakdf_row,"mzmin"]&
             mz<feature_peaks_added[peakdf_row,"mzmax"])
  plot(v$rt, v$int, type="l", 
       main=paste(feature_peaks_added[peakdf_row,"mzmin"], "-", 
                  feature_peaks_added[peakdf_row,"mzmax"]))
  legend("topright", legend=feature_peaks_added[peakdf_row,"qscore"])
}
sorted_feature_peaks <- feature_peaks_added[order(feature_peaks_added$qscore),]
sorted_feature_peaks <- subset(sorted_feature_peaks, qscore>0.1)
for(i in seq_len(nrow(sorted_feature_peaks))){
  plotPeak(rownames(sorted_feature_peaks)[i])
  readline(prompt = "press")
}


# Add MS2 info ----
msms_files <- list.files("mzMLs/MSMS", pattern = "pos", full.names = TRUE)
raw_msmsdata <- lapply(msms_files, grabSingleFileMS2)
nrgs <- as.numeric(gsub(".*neg|.*pos|\\.mzML", "", msms_files))
raw_msmsdata <- lapply(seq_along(raw_msmsdata), function(x){
  cbind(nrg=nrgs[x], raw_msmsdata[[x]])
})
raw_msmsdata <- as.data.table(do.call(rbind, raw_msmsdata))
findMSMSdata <- function(mzr, rtr, ppm=5, rtwindow=10, msmsdata=raw_msmsdata){
  given_msms <- msmsdata[rt%between%c(rtr-rtwindow, rtr+rtwindow)&
                           premz%between%pmppm(mzr, ppm)]
  msms_list <- split(given_msms, given_msms$nrg)
  lapply(msms_list, function(nrg_frags){
    as.data.frame(select(nrg_frags, c("fragmz", "int")))
  })
}
feature_msms <- mapply(FUN = findMSMSdata, mzr=feature_summaries$mz,
                       rtr=feature_summaries$rt, SIMPLIFY = FALSE)
