# Script to run XCMS peakpicking, retention time correction, and grouping
# Called by Control.Rmd
# Requires dev version of XCMS for improved peakpicking given settings

# Read in the raw data ----
# ms_files comes from Control.Rmd
raw_data <- readMSData(files = ms_files, msLevel. = 1, 
                       verbose = TRUE, centroided. = TRUE, mode = "onDisk",
                       pdata = as(metadframe, "AnnotatedDataFrame"))


# Perform peakpicking ----
cwp <- CentWaveParam(ppm = 2.5, peakwidth = c(15, 15), 
                     snthresh = 1, prefilter = c(0, 10000), 
                     integrate = 2, mzCenterFun = "wMean", 
                     mzdiff = 0.001, fitgauss = FALSE, 
                     noise = 5000, firstBaselineCheck = FALSE, 
                     extendLengthMSW = TRUE)
xdata <- suppressMessages(findChromPeaks(raw_data, param = cwp))
saveRDS(xdata, file = paste0(intermediate_folder, "current_xdata.rds"))
message(Sys.time()-start_time)
# 13 minutes


# Assign new quality scores ----
# speedyQscoreCalculator is custom function which should be in functions script
# sourced by speedyQscoreCalculator
xdata <- readRDS(file = paste0(intermediate_folder, "current_xdata.rds"))
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
          file = paste0(intermediate_folder, "all_peaks_w_qscores.csv"))
xdata_cleanpeak <- `chromPeaks<-`(xdata, peakdf_qscored)
saveRDS(xdata_cleanpeak, 
        file = paste0(intermediate_folder, "xdata_cleanpeak.rds"))
message(Sys.time()-start_time)
# 20 minutes



# Other XCMS things (rtcor, group) ----
xdata_cleanpeak <- readRDS(file = paste0(intermediate_folder, "xdata_cleanpeak.rds"))
register(BPPARAM = SnowParam(tasks = length(ms_files), progressbar = TRUE))
obp <- ObiwarpParam(binSize = 0.1, centerSample = 4, 
                    response = 1, distFun = "cor_opt")
xdata_rt <- suppressMessages(adjustRtime(xdata_cleanpeak, param = obp))
plotAdjustedRtime(xdata_rt, col = c("green", "red", "blue", "black", "black")[
  factor(metadframe$depth)])

pdp <- PeakDensityParam(sampleGroups = xdata_rt$depth, 
                        bw = 5, minFraction = 0.5, 
                        binSize = 0.002, minSamples = 2)
xdata_cor <- groupChromPeaks(xdata_rt, param = pdp)

fpp <- FillChromPeaksParam()
xdata_filled <- suppressMessages(fillChromPeaks_wkumler(xdata_cor, param = fpp))

feature_defs <- featureDefinitions(xdata_filled)
raw_peaks <- lapply(seq_len(nrow(feature_defs)), function(i){
  cbind(feature=sprintf("FT%03d", i), 
        peak_id=unlist(feature_defs$peakidx[i]))
}) %>% 
  do.call(what=rbind) %>% 
  as.data.frame(stringsAsFactors=FALSE) %>% 
  mutate(peak_id=as.numeric(peak_id)) %>%
  cbind(chromPeaks(xdata_filled)[.$peak_id, ]) %>%
  mutate(file_name=basename(fileNames(xdata_filled))[sample]) %>%
  arrange(feature, sample)

saveRDS(xdata_filled, file = paste0(intermediate_folder, "current_xdata_filled.rds"))
write.csv(raw_peaks, file = paste0(intermediate_folder, "raw_peaks.csv"), row.names = FALSE)
message(Sys.time()-start_time)
# 30 minutes



