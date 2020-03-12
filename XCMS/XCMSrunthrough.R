
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

ms_files <- "mzMLs" %>%
  list.files(pattern = ".mzML", full.names = TRUE) %>%
  normalizePath() %>%
  `[`(!grepl("Fullneg|Fullpos|QC-KM1906", x = .))



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
xdata_rt <- readRDS(file = "XCMS/temp_data/current_xdata_rt.rds")
pdp <- PeakDensityParam(sampleGroups = xdata_rt$depth, 
                        bw = 2, minFraction = 0.5, 
                        binSize = 0.002)
xdata_cor <- groupChromPeaks(xdata_rt, param = pdp)
saveRDS(xdata_cor, file = "XCMS/temp_data/current_xdata_cor.rds")
print(Sys.time()-start_time)
# 1.8 hours
beep(2)



# Fill peaks ----
xdata_cor <- readRDS(file = "XCMS/temp_data/current_xdata_cor.rds")
xdata_filled <- suppressMessages(fillChromPeaks(xdata_cor, param = FillChromPeaksParam()))
saveRDS(object = xdata_filled, file = "XCMS/temp_data/current_xdata_filled.rds")



# Add new quality scores ----
xdata_filled <- readRDS(file = "XCMS/temp_data/current_xdata_filled.rds")
feature_peaks <- lapply(seq_len(nrow(featureDefinitions(xdata_filled))), function(i){
  cbind(feature=sprintf("FT%03d", i), 
        peak_id=unlist(featureDefinitions(xdata_filled)$peakidx[i]))
}) %>% do.call(what=rbind) %>% 
  as.data.frame(stringsAsFactors=FALSE) %>% 
  mutate(peak_id=as.numeric(peak_id)) %>%
  cbind(chromPeaks(xdata_filled)[.$peak_id, ]) %>%
  mutate(file_name=basename(fileNames(xdata_filled))[sample]) %>%
  select(-c(SNR, peak_cor, qscore))

split_groups <- split(feature_peaks, factor(feature_peaks$file_name))
files_newscores <- bplapply(seq_along(split_groups), 
                            function(x, split_xcms_filedata, ms_files, 
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
split_xcms_filedata=split_groups, 
ms_files=ms_files, 
grabSingleFileData = grabSingleFileData,
xcmsQscoreCalculator = xcmsQscoreCalculator,
qscoreCalculator = qscoreCalculator)
feature_peaks_rescored <- do.call(rbind, files_newscores) %>% 
  as.data.frame() %>% arrange(feature)
saveRDS(feature_peaks_rescored, file = "XCMS/temp_data/feature_peaks_rescored.rds")
print(Sys.time()-start_time)
beep(2)

# Find adducts and isotopes ----
feature_peaks_rescored <- readRDS("XCMS/temp_data/feature_peaks_rescored.rds")
file_peaks <- split(feature_peaks_rescored, factor(feature_peaks_rescored$file_name))
plan(multiprocess)
xdata_adduct_iso <- future_lapply(X = file_peaks, FUN = findAdducts)
xdata_adduct_iso <- do.call(rbind, xdata_adduct_iso)



# Check peak quality (sorry laptop memory)
dt_list <- pblapply(seq_along(ms_files), function(x){
  v <- grabSingleFileData(ms_files[x]) %>%
    cbind(file_name=basename(ms_files[x]), .)
  cor_rt <- unname(split(adjustedRtime(xdata_filled), fromFile(xdata_filled))[[x]])
  v$rt <- factor(v$rt)
  v$cor_rt <- cor_rt[v$rt]
  return(v)
})
dt <- do.call(rbind, dt_list) %>% as.data.table()

while(T){
  ft <- sprintf("FT%03d", i)
  peak_data <- feature_peaks_rescored %>%
    filter(feature==ft) %>%
    summarise(mzmin=min(mzmin), mzmax=max(mzmax), 
              rtmin=min(rtmin), rtmax=max(rtmax),
              m_qscore=mean(qscore), med_qscore=median(qscore))
  eic <- dt[cor_rt%between%c(peak_data$rtmin, peak_data$rtmax)&
              mz%between%c(peak_data$mzmin, peak_data$mzmax)] %>%
    mutate(file_name=as.character(file_name))
  mdframe <- metadata@data %>% 
    mutate(filenames=paste0("190715_", as.character(filenames)))
  eic_df <- left_join(eic, mdframe, by=c("file_name"="filenames"))
  gp <- ggplot(eic_df) + 
    geom_line(aes(x=cor_rt, y=int, group=file_name, color=depth)) +
    scale_color_manual(values = c(Blank="red", Pooled="black", `25m`="blue", 
                                  DCM="green", Std="gray")) +
    ggtitle(paste("Mean quality:", round(peak_data$m_qscore), 
                  "\nMed quality:", round(peak_data$med_qscore)))
    
  print(gp)
  readline(prompt = "Press Enter")
  i <- i+1
}

peak_data <- feature_peaks_rescored %>%
  group_by(feature) %>%
  summarise(m_qscore=mean(qscore), med_qscore=median(qscore))
summary(peak_data)
