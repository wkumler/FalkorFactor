
# Setup things ----
library(xcms)
library(tidyverse)
library(data.table)
library(beepr)
library(httr)
library(future.apply)
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
findAdducts <- function(x){
  file_data <- fileNames(xdata_cor)[unique(x$sample)] %>%
    grabSingleFileData() %>%
    as.data.table()
  file_data$rt <- adjustedRtime(xdata_cor) %>%
    split(fromFile(xdata_cor)) %>%
    `[[`(unique(x$sample)) %>%
    `[`(factor(file_data$rt))
  outlist <- list()
  for(i in seq_len(nrow(x))){
    peak_row_data <- x[i, ]
    init_eic <- file_data[mz %between% c(peak_row_data$mzmin, peak_row_data$mzmax)&
                            rt %between% c(peak_row_data$rtmin, peak_row_data$rtmax)]
    mp1 <- checkCor(mass = peak_row_data$mz+1.003355, rtmin = peak_row_data$rtmin, 
                    rtmax = peak_row_data$rtmax, init_eic = init_eic, file_data=file_data)
    mp2 <- checkCor(mass = peak_row_data$mz+1.003355*2, rtmin = peak_row_data$rtmin, 
                    rtmax = peak_row_data$rtmax, init_eic = init_eic, file_data=file_data)
    m_na <- checkCor(mass = peak_row_data$mz-1.007276+22.98922, init_eic = init_eic,
                     rtmin = peak_row_data$rtmin, rtmax = peak_row_data$rtmax, file_data=file_data)
    m_h <- checkCor(mass = peak_row_data$mz-22.98922+1.007276, init_eic = init_eic,
                    rtmin = peak_row_data$rtmin, rtmax = peak_row_data$rtmax, file_data=file_data)
    outlist[[i]] <- data.frame(M_1_area = mp1[1], M_1_q = mp1[2], M_1_simil = mp1[3],
                               M_2_area = mp2[1], M_2_q = mp2[2], M_2_simil = mp2[3],
                               M_Na_area = m_na[1], M_Na_q = m_na[2], M_Na_simil = m_na[3],
                               M_H_area = m_h[1], M_H_q = m_h[2], M_H_simil = m_h[3])
  }
  return(do.call(rbind, outlist))
}
checkCor <- function(mass, rtmin, rtmax, init_eic, file_data){
  given_eic <- file_data[mz %between% pmppm(mass, ppm = 5) & rt %between% c(rtmin, rtmax)]
  if(nrow(given_eic)<3){
    return(c(0L, 0L, 0L))
  }
  peak_qscore <- qscoreCalculator(given_eic)$qscore
  merged_eic <- merge(init_eic, given_eic, by="rt")
  if(nrow(given_eic)<3){
    return(c(0L, 0L, 0L))
  }
  peak_match <- cor(merged_eic$int.x, merged_eic$int.y)
  peak_area <- trapz(given_eic$rt, given_eic$int)
  return(c(peak_area, peak_qscore, peak_match))
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

register(BPPARAM = SnowParam(tasks = length(ms_files), progressbar = TRUE))

raw_data <- readMSData(files = ms_files, pdata = metadata, mode = "onDisk")
saveRDS(raw_data, file = "XCMS/temp_data/current_raw_data.rds")
print(Sys.time()-start_time)
#2.3 minutes
beep(2)



# Perform peakpicking ----
raw_data <- readRDS(file = "XCMS/temp_data/current_raw_data.rds")
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
peakdf_qscored <- as.matrix(read.csv(file = "XCMS/temp_data/peakdf_qscored.csv"))
threshold <- 20
cleandf_qscored <- peakdf_qscored[peakdf_qscored[, "qscore"]>threshold,]
xdata_cleanpeak <- `chromPeaks<-`(xdata, cleandf_qscored)
# sample_df <- peakdf_qscored[sample(1:nrow(peakdf_qscored), size = 10000),]
# plot(log10(sample_df[,"sn"]), log10(sample_df[,"qscore"]))

# Adjust retention time and compare ----
obp <- ObiwarpParam(binSize = 0.01, centerSample = 4, response = 1, distFun = "cor_opt")
xdata_rt <- suppressMessages(adjustRtime(xdata_cleanpeak, param = obp))
plotAdjustedRtime(xdata_rt)
saveRDS(xdata_rt, file = "XCMS/temp_data/current_xdata_rt.rds")
print(Sys.time()-start_time)
# 1.7 hours
beep(2)



# Correspondence ----
xdata_rt <- readRDS("XCMS/temp_data/current_xdata_rt.rds")
pdp <- PeakDensityParam(sampleGroups = xdata_rt$depth, 
                        bw = 5, minFraction = 0.5, 
                        binSize = 0.002)
xdata_cor <- groupChromPeaks(xdata_rt, param = pdp)
saveRDS(xdata_cor, file = "XCMS/temp_data/current_xdata_cor.rds")
print(Sys.time()-start_time)
# 1.8 hours
beep(2)



# Fill peaks and add new quality scores ----
xdata_cor <- readRDS(file = "XCMS/temp_data/current_xdata_cor.rds")
xdata_filled <- suppressMessages(fillChromPeaks(xdata_cor, param = FillChromPeaksParam()))
peak_df <- as.data.frame(chromPeaks(xdata_filled))
peak_list <- split(peak_df, is.na(peak_df$qscore))
unscored_rows <- peak_list[["TRUE"]] %>%
  select(-c("SNR", "peak_cor", "qscore"))
split_unscored <- split(unscored_rows, unscored_rows$sample)
fileids <- unique(peak_df$sample)
files_newscores <- bplapply(fileids, function(x, split_xcms_filedata, ms_files, 
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
split_xcms_filedata=split_unscored, 
ms_files=ms_files, 
grabSingleFileData = grabSingleFileData,
xcmsQscoreCalculator = xcmsQscoreCalculator,
qscoreCalculator = qscoreCalculator)
peak_list[["TRUE"]] <- do.call(rbind, files_newscores)
chromPeaks(xdata_filled) <- `rownames<-`(as.matrix(do.call(rbind, peak_list)), rownames(peak_df))
pdp <- PeakDensityParam(sampleGroups = xdata_filled$depth, 
                        bw = 5, minFraction = 0.5, 
                        binSize = 0.002)
xdata_filled <- groupChromPeaks(xdata_filled, param = pdp)
saveRDS(object = xdata_filled, file = "XCMS/temp_data/current_xdata_filled.rds")


# Find adducts and isotopes ----
xdata_filled <- readRDS(file = "XCMS/temp_data/current_xdata_filled.rds")

peak_df <- as.data.frame(chromPeaks(xdata_filled))
group_df <- featureDefinitions(xdata_filled) %>%
  cbind(featureValues(xdata_filled)) %>%
  as.data.frame() %>%
  mutate(feature_name = rownames(.)) %>%
  pivot_longer(cols = starts_with("X190715"), names_to = "file_name", values_to = "peak_area")

unduped_peak_idxs <- lapply(seq_len(nrow(featureDefinitions(xdata_filled))), function(i){
  peak_idxs <- unlist(featureDefinitions(xdata_filled)[i,"peakidx"])
  if(length(peak_idxs==1)){
    return(peak_idxs)
  }
  subs <- as.data.frame(chromPeaks(xdata_filled)[peak_idxs, ]) %>% 
    mutate(peakidx=peak_idxs) %>%
    group_by(sample) %>% 
    slice(which.max(qscore)) %>%
    pull(peakidx)
  return(subs)
})


# Conundrum: 
group_df %>% group_by(feature_name) %>% summarize(sum(!is.na(peak_area))) %>% head()
head(sapply(featureDefinitions(xdata_filled)$peakidx, length))



file_peaks <- split(peak_df, peak_df$sample)
plan(multiprocess)
xdata_adduct_iso <- future_lapply(X = file_peaks, FUN = findAdducts)
xdata_adduct_iso <- do.call(rbind, xdata_adduct_iso)
xdata_adduct_iso <- cbind(peak_df, xdata_adduct_iso)




# Investigate some peaks! ----
xdata_cor <- readRDS(file = "XCMS/temp_data/current_xdata_cor.rds")
peak_idx_df <- chromPeaks(xdata_cor)
cor_peaklist <- featureDefinitions(xdata_cor) %>%
  cbind(featureValues(xdata_cor)) %>%
  as.data.frame() %>%
  filter(rtmax>60&rtmin<1100) %>%
  mutate(coll_min_rt=sapply(peakidx, function(x)min(peak_idx_df[x, "rtmin"]))) %>%
  mutate(coll_max_rt=sapply(peakidx, function(x)max(peak_idx_df[x, "rtmax"]))) %>%
  mutate(coll_min_mz=sapply(peakidx, function(x)min(peak_idx_df[x, "mzmin"]))) %>%
  mutate(coll_max_mz=sapply(peakidx, function(x)max(peak_idx_df[x, "mzmax"])))


rt_cors <- split(unname(adjustedRtime(xdata_cor)), fromFile(xdata_cor))
names(rt_cors) <- basename(fileNames(xdata_cor))

blank_data <- grabSingleFileData(filename = fileNames(xdata_cor)[1])
blank_data$cor_rt <- rt_cors[[1]][factor(blank_data$rt)]
blank_data <- as.data.table(blank_data)

blank_peaks <- cor_peaklist %>%
  filter(cor_peaklist$Blank==1&npeaks==1) %>%
  select(Blank, starts_with("coll"), X190715_Blk_KM1906U14.Blk_C.mzML)
for(i in seq_len(nrow(blank_peaks))){
  peak_row_data <- blank_peaks[i,]
  eic <- blank_data[cor_rt %between% c(peak_row_data$coll_min_rt, peak_row_data$coll_max_rt)&
                    mz %between% c(peak_row_data$coll_min_mz, peak_row_data$coll_max_mz)]
  plot(eic$cor_rt, eic$int, type = "l")
  readline(prompt = "Press Enter")
}
masses <- (blank_peaks$coll_min_mz+blank_peaks$coll_max_mz)/2-1.007276
pubchemresults <- lapply(masses, searchPubchem, ppm = 5)
blank_peaks$pubchem <- pubchemresults
head(blank_peaks)




pool_data <- grabSingleFileData(filename = fileNames(xdata_cor)[3])
pool_data$cor_rt <- rt_cors[[3]][factor(pool_data$rt)]
pool_data <- as.data.table(pool_data)

betaine_na <- 117.078979+22.989222
betaine_h <- 117.078979+1.007276
betaine_d <- betaine_h+1.003355
betaine_d2 <- betaine_h+1.003355*2

betaine_na_eic <- pool_data[mz %between% pmppm(betaine_na, ppm=5)]
betaine_h_eic <- pool_data[mz %between% pmppm(betaine_h, ppm=5)]
betaine_d_eic <- pool_data[mz %between% pmppm(betaine_d, ppm=5)]
betaine_d2_eic <- pool_data[mz %between% pmppm(betaine_d2, ppm=5)]

plot(betaine_na_eic$rt, betaine_na_eic$int, type="l", xlim=c(300, 400))
plot(betaine_h_eic$rt, betaine_h_eic$int, type="l", xlim=c(300, 400))
plot(betaine_d_eic$rt, betaine_d_eic$int, type="l", xlim=c(300, 400))
plot(betaine_d2_eic$rt, betaine_d2_eic$int, type="l", xlim=c(300, 400))

v <- featureDefinitions(xdata_filled)
bet_found <- v[betaine_h>v$mzmin&v$h<v$mzmax]