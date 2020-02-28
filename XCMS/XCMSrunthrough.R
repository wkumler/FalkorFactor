
# Setup things ----
library(xcms)
library(dplyr)
library(RSQLite)
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
  if(nrow(eic)<3){
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
  return(SNR*peak_cor^4*log10(max(eic$int)))
}
xcmsQscoreCalculator <- function(df_row, xcms_peakdf, database){
  #Extract the relevant EIC
  peak_row_data <- xcms_peakdf[df_row, ]
  eic <- dbGetQuery(conn=database, paste(
  "SELECT * FROM MS1_pos",
  paste0("WHERE filename = '", basename(ms_files[peak_row_data$sample]), "'"),
  "AND mz >=", peak_row_data$mzmin,
  "AND mz <=", peak_row_data$mzmax,
  "AND rt >=", peak_row_data$rtmin,
  "AND rt <=", peak_row_data$rtmax
  ))
  dbClearResult(eic)
  return(qscoreCalculator(eic))
}



# Load MS data ----
ms_files <- "../mzMLs" %>%
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



# Create SQLite database for rapid manual retrieval ----
if(file.exists("temp_data/falkor.db"))file.remove("temp_data/falkor.db")
if(!file.exists("temp_data/falkor.db"))file.create("temp_data/falkor.db")

falkor_db <- dbConnect(drv = RSQLite::SQLite(), "temp_data/falkor.db")
rs <- dbSendQuery(conn=falkor_db, "CREATE TABLE MS1_pos (
            filename TEXT, mz NUMERIC, rt NUMERIC, int NUMERIC)")
dbClearResult(rs)
for(i in ms_files){
  filedata <- grabSingleFileData(i)
  filedata$filename <- basename(i)
  print(head(filedata))
  dbWriteTable(conn=falkor_db, name="MS1_pos", filedata, append=T, row.names=F)
}
rs <- dbSendQuery(conn = falkor_db, "CREATE INDEX mz_pos_lookup ON MS1_pos(mz)")
dbClearResult(rs)
rs <- dbSendQuery(conn = falkor_db, "CREATE INDEX rt_pos_lookup ON MS1_pos(rt)")
dbClearResult(rs)
rs <- dbSendQuery(conn = falkor_db, "CREATE INDEX filename_pos_lookup ON MS1_pos(filename)")
dbClearResult(rs)
dbDisconnect(falkor_db)
print(Sys.time()-start_time)


# Perform peakpicking ----
raw_data <- readRDS(file = "temp_data/current_raw_data.rds")
cwp <- CentWaveParam(ppm = 5, peakwidth = c(20, 80), 
                     snthresh = 0, prefilter = c(0, 0), 
                     integrate = 1, mzCenterFun = "wMean", 
                     mzdiff = 0.0001, fitgauss = FALSE, 
                     noise = 0, firstBaselineCheck = FALSE)
xdata <- findChromPeaks(raw_data, param = cwp)
print(xdata)
saveRDS(xdata, file = "temp_data/current_xdata.rds")
print(Sys.time()-start_time)
# 34.2 minutes
beep(2)



# Re-assign quality scores to confirm good peaks FIX THIS ----
xdata <- readRDS(file = "temp_data/current_xdata.rds")
xcms_peakdf <- as.data.frame(chromPeaks(xdata))

# Consider instead, splitting by file and loading in one file's raw data at a time
# via an lapply() that then calls the qscore. Nicely parallelizable!
falkor_db <- dbConnect(drv = RSQLite::SQLite(), "temp_data/falkor.db")
qscores <- sapply(seq_len(nrow(xcms_peakdf)), xcmsQscoreCalculator,
                  xcms_peakdf=xcms_peakdf, database=falkor_db)
dbDisconnect(falkor_db)


# Get peaklist and poke around ----
xdata <- readRDS(file = "temp_data/current_xdata.rds")
xcms_peaks <- as.data.frame(chromPeaks(xdata))
sum(xcms_peaks$sn==round(xcms_peaks$maxo)-1) # Number of weird S/N estimates
xcms_peaks %>% group_by(sample) %>% summarise(peak_count=n()) # Number of peaks in each

arsb_rtr <- c(300, 420)
arsb_mz <- 177.997499+1.007276
arsb_mzr <- pmppm(arsb_mz, ppm = 5)
chr_raw <- chromatogram(xdata, rt = arsb_rtr, mz = arsb_mzr)
plot(chr_raw)
plot(chr_raw, ylim=c(0, 4000000))



# Adjust retention time and compare ----
obp <- ObiwarpParam(binSize = 0.01, centerSample = 4, response = 1, distFun = "cor_opt")
xdata_rt <- adjustRtime(xdata, param = obp)
plotAdjustedRtime(xdata_rt)
saveRDS(xdata_rt, file = "temp_data/current_xdata_rt.rds")
print(Sys.time()-start_time)
# 50 minutes
beep(2)

chr_adj <- chromatogram(xdata_rt, rt = arsb_rtr, mz = arsb_mzr)
plot(chr_adj)
plot(chr_adj, ylim=c(0, 4000000))
arsb_adj <- do.call(rbind, lapply(chr_adj, function(x)x@chromPeaks))
cbind(sapply(X = chr_adj, function(x)nrow(x@chromPeaks)), basename(fileNames(xdata_rt)))
boxplot(arsb_adj[,"rt"])



# Correspondence ----
xdata_rt <- readRDS("temp_data/current_xdata_rt.rds")
pdp <- PeakDensityParam(sampleGroups = xdata_rt$sample_group, 
                        bw = 20, minFraction = 0.5, 
                        binSize = 0.002)
xdata_cor <- groupChromPeaks(xdata_rt, param = pdp)
print(Sys.time()-start_time)
beep(2)


# Group isotopes and adducts ----

