# Script to run XCMS peakpicking, retention time correction, and grouping
# Called by Control.Rmd
# Requires dev version of XCMS for improved peakpicking given settings

# ms_files comes from Control.Rmd
# qscore_threshold comes from Control.Rmd


# Read in the raw data ----
start_time <- Sys.time()
raw_data <- readMSData(files = normalizePath(paste("mzMLs", polarity, ms_files, sep = "/")), 
                       msLevel. = 1, 
                       verbose = TRUE, centroided. = TRUE, mode = "onDisk",
                       pdata = as(metadframe, "AnnotatedDataFrame"))
message("Time to read in files: ", round(Sys.time()-start_time, digits = 2), " min")

# Perform peakpicking ----
start_time <- Sys.time()
cwp <- CentWaveParam(ppm = 2.5, peakwidth = c(15, 15), 
                     snthresh = 1, prefilter = c(0, 10000), 
                     integrate = 2, mzCenterFun = "wMean", 
                     mzdiff = 0.001, fitgauss = FALSE, 
                     noise = 5000, firstBaselineCheck = FALSE, 
                     extendLengthMSW = TRUE)
xdata <- suppressMessages(findChromPeaks(raw_data, param = cwp))
message("Time to perform peakpicking: ", 
        round(Sys.time()-start_time, digits = 2), " min")
# 13 minutes


# Assign new quality scores ----
# speedyQscoreCalculator is custom function which should be in functions script
# sourced by speedyQscoreCalculator
start_time <- Sys.time()
xcms_peakdf <- chromPeaks(xdata) %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  mutate(filename=fileNames(xdata)[.[["sample"]]]) %>%
  arrange(mz)
split_xcms_filedata <- split(xcms_peakdf, xcms_peakdf$filename)
files_qscores <- bplapply(X = split_xcms_filedata, FUN = speedyQscoreCalculator,
                          grabSingleFileData, qscoreCalculator)

peakdf_qscored <- files_qscores %>%
  do.call(what = rbind) %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  filter(sn>qscore_threshold) %>%
  select(-filename) %>%
  arrange(sample, rtmin, rtmax) %>%
  as.matrix()
xdata_cleanpeak <- `chromPeaks<-`(xdata, peakdf_qscored)
message("Time to assign quality scores: ", 
        round(Sys.time()-start_time, digits = 2), " min")
# 20 minutes



# Other XCMS things (rtcor, group) ----
start_time <- Sys.time()
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

message("Time for other XCMS things: ", 
        round(Sys.time()-start_time, digits = 2), " min")
# 30 minutes
