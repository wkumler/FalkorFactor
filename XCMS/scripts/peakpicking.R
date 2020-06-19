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
                        checkPeakCor, pmppm, trapz, polarity){
  file_peaks <- file_peaks[order(file_peaks$rtmax),]
  file_path <- paste("mzMLs/", polarity, unique(file_peaks$file_name), sep = "/")
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
  given_eic <- file_data_table[mz%between%pmppm(mass, ppm = 5) & 
                                 rt%between%c(rtmin, rtmax)]
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
                          checkPeakCor, pmppm, trapz, polarity){
  file_peaks <- file_peaks[order(file_peaks$rtmax),]
  file_path <- paste("mzMLs", polarity, unique(file_peaks$file_name), sep = "/")
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
fillChromPeaks_wkumler <- function(object, param){
  msLevel <- 1L
  BPPARAM <- bpparam()
  
  if (length(msLevel) != 1) 
    stop("Can only perform peak filling for one MS level at a time")
  if (!hasFeatures(object, msLevel = msLevel)) 
    stop("No feature definitions for MS level ", 
         msLevel, " present. Please run 'groupChromPeaks' first.")
  if (xcms:::.hasFilledPeaks(object)) 
    message("Filled peaks already present, adding still missing", 
            " peaks.")
  if (xcms:::hasChromPeaks(object) & !xcms:::.has_chrom_peak_data(object)) 
    object <- updateObject(object)
  startDate <- date()
  expandMz <- expandMz(param)
  expandRt <- expandRt(param)
  fixedMz <- fixedMz(param)
  fixedRt <- fixedRt(param)
  ppm <- ppm(param)
  message("Defining peak areas for filling-in .", 
          appendLF = FALSE)
  fdef <- featureDefinitions(object, msLevel = msLevel)
  aggFunLow <- median
  aggFunHigh <- median
  tmp_pks <- chromPeaks(object)[, c("rtmin", "rtmax", 
                                    "mzmin", "mzmax")]
  pkArea <- do.call(rbind, lapply(fdef$peakidx, function(z) {
    pa <- c(aggFunLow(tmp_pks[z, 1]), aggFunHigh(tmp_pks[z, 2]), 
            aggFunLow(tmp_pks[z, 3]), aggFunHigh(tmp_pks[z, 4]))
    if (ppm != 0) {
      mzmean <- mean(pa[3:4])
      tittle <- mzmean * (ppm/2)/1e+06
      if ((pa[4] - pa[3]) < (tittle * 2)) {
        pa[3] <- mzmean - tittle
        pa[4] <- mzmean + tittle
      }
    }
    if (expandRt != 0) {
      diffRt <- (pa[2] - pa[1]) * expandRt/2
      pa[1] <- pa[1] - diffRt
      pa[2] <- pa[2] + diffRt
    }
    if (expandMz != 0) {
      diffMz <- (pa[4] - pa[3]) * expandMz/2
      pa[3] <- pa[3] - diffMz
      pa[4] <- pa[4] + diffMz
    }
    if (fixedMz != 0) {
      pa[3] <- pa[3] - fixedMz
      pa[4] <- pa[4] + fixedMz
    }
    if (fixedRt != 0) {
      pa[1] <- pa[1] - fixedRt
      pa[2] <- pa[2] + fixedRt
    }
    pa
  }))
  rm(tmp_pks)
  message(".", appendLF = FALSE)
  colnames(pkArea) <- c("rtmin", "rtmax", "mzmin", "mzmax")
  pkArea <- cbind(group_idx = 1:nrow(pkArea), pkArea, mzmed = as.numeric(fdef$mzmed))
  pkGrpVal <- featureValues(object, value = "index", msLevel = msLevel)
  message(".", appendLF = FALSE)
  if (!any(is.na(rowSums(pkGrpVal)))) {
    message("No missing peaks present.")
    return(object)
  }
  message(".", appendLF = FALSE)
  objectL <- vector("list", length(fileNames(object)))
  pkAreaL <- objectL
  req_fcol <- requiredFvarLabels("OnDiskMSnExp")
  min_fdata <- fData(object)[, req_fcol]
  rt_range <- range(pkArea[, c("rtmin", "rtmax")])
  if (hasAdjustedRtime(object)) 
    min_fdata$retentionTime <- adjustedRtime(object)
  for (i in 1:length(fileNames(object))) {
    fd <- min_fdata[min_fdata$fileIdx == i, ]
    fd$fileIdx <- 1L
    objectL[[i]] <- new("OnDiskMSnExp", 
                        processingData = new("MSnProcess", 
                                             files = fileNames(object)[i]), 
                        featureData = new("AnnotatedDataFrame", fd), 
                        phenoData = new("NAnnotatedDataFrame", 
                                        data.frame(sampleNames = "1")), 
                        experimentData = new("MIAPE", 
                                             instrumentManufacturer = "a", 
                                             instrumentModel = "a", 
                                             ionSource = "a", 
                                             analyser = "a", 
                                             detectorType = "a"))
    pkAreaL[[i]] <- pkArea[is.na(pkGrpVal[, i]), , drop = FALSE]
  }
  rm(pkGrpVal)
  rm(pkArea)
  rm(min_fdata)
  message(" OK\nStart integrating peak areas from original files")
  ph <- processHistory(object, type = xcms:::.PROCSTEP.PEAK.DETECTION)
  findPeakMethod <- "unknown"
  mzCenterFun <- "wMean"
  if (length(ph)) {
    if (is(ph[[1]], "XProcessHistory")) {
      prm <- ph[[1]]@param
      findPeakMethod <- .param2string(prm)
      if (.hasSlot(prm, "mzCenterFun")) 
        mzCenterFun <- prm@mzCenterFun
    }
  }
  cp_colnames <- colnames(chromPeaks(object))
  mzCenterFun <- paste("mzCenter", gsub(mzCenterFun, 
                                        pattern = "mzCenter.", replacement = "", 
                                        fixed = TRUE), sep = ".")
  if (findPeakMethod == "MSW") {
    rts <- rtime(object, bySample = TRUE)
    if (any(lengths(rts) > 1)) 
      stop("The data is supposed to be direct injection data, ", 
           "but I got files with more than one spectrum/", 
           "retention time!")
    res <- bpmapply(FUN = xcms:::.getMSWPeakData, objectL, pkAreaL, 
                    as.list(1:length(objectL)), MoreArgs = list(cn = cp_colnames), 
                    BPPARAM = BPPARAM, SIMPLIFY = FALSE)
  } else if (findPeakMethod == "matchedFilter") {
    res <- bpmapply(FUN = xcms:::.getChromPeakData_matchedFilter, 
                    objectL, pkAreaL, as.list(1:length(objectL)), 
                    MoreArgs = list(cn = cp_colnames, param = prm, 
                                    msLevel = msLevel), 
                    BPPARAM = BPPARAM, 
                    SIMPLIFY = FALSE)
  } else {
    res <- bpmapply(FUN = xcms:::.getChromPeakData, objectL, 
                    pkAreaL, as.list(1:length(objectL)), 
                    MoreArgs = list(cn = cp_colnames, mzCenterFun = mzCenterFun, 
                                    msLevel = msLevel), 
                    BPPARAM = BPPARAM, SIMPLIFY = FALSE)
  }
  rm(objectL)
  res <- do.call(rbind, res)
  res <- cbind(res, group_idx = unlist(lapply(pkAreaL, function(z){
    z[, "group_idx"]}), use.names = FALSE))
  res <- res[!is.na(res[, "into"]), , drop = FALSE]
  if (nrow(res) == 0) {
    warning("Could not integrate any signal for the missing ", 
            "peaks! Consider increasing 'expandMz' and 'expandRt'.")
    return(object)
  }
  rm(pkAreaL)
  gc()
  newFd <- new("MsFeatureData")
  newFd@.xData <- xcms:::.copy_env(object@msFeatureData)
  object@msFeatureData <- new("MsFeatureData")
  incr <- nrow(chromPeaks(newFd))
  for (i in unique(res[, "group_idx"])) {
    fdef$peakidx[[i]] <- c(fdef$peakidx[[i]], (which(res[, "group_idx"] == i) + incr))
  }
  fdef <- rbind(fdef, featureDefinitions(newFd)[featureDefinitions(newFd)$ms_level != 
                                                  msLevel, , drop = FALSE])
  if (!any(colnames(fdef) == "ms_level")) {
    fdef$ms_level <- 1L
  } else {fdef <- fdef[order(fdef$ms_level), ]}
  maxId <- max(as.numeric(gsub("^.*\\.", "", rownames(chromPeaks(newFd)))))
  if (maxId < 1)
    stop("chromPeaks matrix lacks rownames; please update ", 
         "'object' with the 'updateObject' function.")
  toId <- maxId + nrow(res)
  rownames(res) <- sprintf(paste0("CP", "%0", ceiling(log10(toId + 1L)), "d"), 
                           (maxId + 1L):toId)
  chromPeaks(newFd) <- rbind(chromPeaks(newFd), res[, -ncol(res)])
  cpd <- chromPeakData(newFd)[rep(1L, nrow(res)), , drop = FALSE]
  cpd[, ] <- NA
  cpd$ms_level <- as.integer(msLevel)
  cpd$is_filled <- TRUE
  if (!any(colnames(chromPeakData(newFd)) == "is_filled")) 
    chromPeakData(newFd)$is_filled <- FALSE
  chromPeakData(newFd) <- rbind(chromPeakData(newFd), cpd)
  rownames(chromPeakData(newFd)) <- rownames(chromPeaks(newFd))
  featureDefinitions(newFd) <- fdef
  lockEnvironment(newFd, bindings = TRUE)
  object@msFeatureData <- newFd
  ph <- xcms:::XProcessHistory(param = param, date. = startDate, 
                               type. = xcms:::.PROCSTEP.PEAK.FILLING, 
                               fileIndex = 1:length(fileNames(object)), 
                               msLevel = msLevel)
  object <- xcms:::addProcessHistory(object, ph)
  object
}



# Setup things ----
start_time <- Sys.time()
library(tidyverse)
library(data.table)
library(pbapply)
library(xcms)

polarity <- "pos"
#polarity <- "neg"
pretty_folder <- paste0("XCMS/", polarity, "_pretty/")
intermediate_folder <- paste0("XCMS/", polarity, "_intermediate/")

ms_files <- "mzMLs" %>%
  paste0("/", polarity) %>%
  list.files(pattern = ".mzML", full.names = TRUE) %>%
  normalizePath() %>%
  `[`(!grepl("Fullneg|Fullpos|QC-KM1906", x = .))

metadframe <- data.frame(
  fileid=basename(ms_files),
  depth=regmatches(regexpr(pattern = "Std|Poo|Blk|DCM|25m", ms_files), x = ms_files),
  station=regmatches(regexpr(pattern = "Std|Poo|Blk|S62|S64|S77|S80", ms_files), 
                     x = ms_files)
)
station_spindirs <- c(Blk="Blk", S62="Cyclone", S64="Cyclone", 
                      S77="Anticyclone", S80="Anticyclone",
                      Poo="Poo", Std="Std")
metadframe$spindir=station_spindirs[metadframe$station]
station_times <- c(Blk="Blk", S62="Morning", S64="Afternoon", 
                   S77="Morning", S80="Afternoon",
                   Poo="Poo", Std="Std")
metadframe$time=station_times[metadframe$station]
write.csv(x = metadframe, row.names = FALSE,
          file = paste0(pretty_folder, "metadata.csv"))

register(BPPARAM = SnowParam(tasks = length(ms_files), progressbar = TRUE))


# Peakpicking ----
raw_data <- readMSData(files = ms_files, msLevel. = 1, 
                       verbose = TRUE, centroided. = TRUE, mode = "onDisk",
                       pdata = as(metadframe, "AnnotatedDataFrame"))
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



# Find isotopes and adducts ----
xdata_filled <- readRDS(file = paste0(intermediate_folder, "current_xdata_filled.rds"))
raw_peaks <- read.csv(file = paste0(intermediate_folder, "raw_peaks.csv"))

is_peak_iso <- 
  bplapply(split(raw_peaks, raw_peaks$file_name),
           FUN = isIsoAdduct, xdata=xdata_filled,
           grabSingleFileData=grabSingleFileData, checkPeakCor=checkPeakCor, 
           pmppm=pmppm, trapz=trapz, polarity=polarity) %>%
  do.call(what = rbind) %>% as.data.frame()
write.csv(is_peak_iso, file = paste0(intermediate_folder, "is_peak_iso.csv"), row.names = FALSE)
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
}) %>% 
  do.call(what=rbind) %>% 
  `[<-`(is.na(.), 0) %>% 
  as.data.frame(stringsAsFactors=FALSE) %>%
  mutate(feature=unique(is_peak_iso$feature)) %>%
  select(feature, everything()) %>%
  arrange(feature)

likely_addisos <- peakareamatch$feature[
  which(rowSums(peakareamatch[,names(peakareamatch)!="feature"]>0.95&
                  peakshapematch[,names(peakshapematch)!="feature"]>0.9)>=1)
  ]

addiso_features <- raw_peaks %>%
  group_by(feature) %>%
  summarise(mzmed=median(mz), rtmed=median(rt), avginto=mean(into, na.rm=TRUE)) %>%
  filter(feature%in%likely_addisos)
write.csv(addiso_features, file = paste0(pretty_folder, "addiso_features.csv"), 
          row.names = FALSE)

message(Sys.time()-start_time)
#40 minutes


# Calculate isotopes and adducts for remaining peaks ----
raw_peaks <- read.csv(paste0(intermediate_folder, "raw_peaks.csv"))
addiso_features <- read.csv(paste0(pretty_folder, "addiso_features.csv"))

# For each peak, look for data at +/- each adduct/isotope m/z 
# Also calculate cor while the raw data is being accessed anyway
complete_peaks <- raw_peaks %>%
  filter(!feature%in%addiso_features$feature) %>%
  split(.$file_name) %>%
  bplapply(FUN = findIsoAdduct, xdata=xdata_filled,
           grabSingleFileData=grabSingleFileData, checkPeakCor=checkPeakCor, 
           pmppm=pmppm, trapz=trapz, polarity=polarity) %>%
  do.call(what = rbind) %>% as.data.frame() %>% 
  `rownames<-`(NULL) %>% arrange(feature)

# Calculate median cor for each FEATURE from the various peak cors
peak_cors <- complete_peaks %>%
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
peak_slope_R2 <- lapply(unique(complete_peaks$feature), function(i){
  feature_areas <- complete_peaks[complete_peaks$feature==i,]
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
  mutate(feature=unique(complete_peaks$feature)) %>%
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
complete_features <- complete_peaks %>% 
  group_by(feature) %>%
  summarize(mzmed=median(mz), rtmed=median(rt), avgarea=mean(M_area)) %>%
  left_join(peak_cors, by="feature") %>%
  left_join(peak_R2s, by=c("feature", "addiso")) %>%
  left_join(peak_slopes, by=c("feature", "addiso")) %>%
  mutate(rel_int=ifelse(cor>0.8&R2>0.9, round(slope*avgarea), 0)) %>%
  select(-c("cor", "R2", "slope")) %>%
  pivot_wider(names_from = addiso, values_from = rel_int)

write.csv(x = complete_peaks, 
          file = paste0(intermediate_folder, "complete_peaks.csv"),
          row.names = FALSE)
write.csv(x = complete_features, 
          file = paste0(intermediate_folder, "complete_features.csv"),
          row.names = FALSE)



# Normalize to the best internal standard ----
xdata_filled <- readRDS(paste0(intermediate_folder, "current_xdata_filled.rds"))
raw_peaks <- read.csv(paste0(intermediate_folder, "raw_peaks.csv"))
addiso_features <- read.csv(paste0(pretty_folder, "addiso_features.csv"))
is_peak_iso <- read.csv(paste0(intermediate_folder, "is_peak_iso.csv"))
complete_peaks <- read.csv(paste0(intermediate_folder, "complete_peaks.csv"))
bionorm_values <- "XCMS/Sample.Key.Falkor.Manual.csv" %>%
  read.csv() %>%
  select(file_name="ï..Sample.Name", norm_vol="Bio.Normalization")
cut.off <- 0.4 #Necessary improvement for "acceptable"
cut.off2 <- 0.1 #If RSD already below, skip B-MIS

# Grab the internal standards and clean up a little
internal_stans <- read.csv("XCMS/falkor_stans.csv") %>%
  filter(Compound.Type=="Internal Standard") %>%
  filter(Fraction1==paste0("HILIC", paste0(toupper(substring(polarity, 1, 1)), 
                                           substring(polarity, 2)))) %>%
  mutate(m.z=as.numeric(m.z)) %>%
  mutate(lower_mz_bound=lapply(m.z, pmppm, ppm=5) %>% sapply(`[`, 1)) %>%
  mutate(upper_mz_bound=lapply(m.z, pmppm, ppm=5) %>% sapply(`[`, 2)) %>%
  mutate(RT_sec=RT..min.*60) %>%
  select(Compound.Name, Emperical.Formula, RT_sec, 
         m.z, lower_mz_bound, upper_mz_bound)

# For each IS, look in the picked peaks and see if one matches mz & rt
found_stans <- internal_stans %>%
  split(seq_len(nrow(.))) %>%
  lapply(function(i){
    suppressMessages(
      raw_peaks %>%
        group_by(feature) %>%
        summarize(mzmed=median(mz), rtmed=median(rt)) %>%
        filter(mzmed%between%c(i$lower_mz_bound,i$upper_mz_bound)) %>%
        mutate(stan=i$Compound.Name) %>%
        mutate(ppm_diff=(abs(i$m.z-.$mzmed)/.$mzmed)*1000000) %>%
        mutate(rt_diff=i$RT_sec-.$rtmed) %>%
        select(feature, stan, mzmed, ppm_diff, rtmed, rt_diff)
    )
  }) %>% do.call(what="rbind") %>% as.data.frame()

# Remove all stans for which more than one peak was found
found_stans <- found_stans[
  !found_stans$stan%in%found_stans$stan[duplicated(found_stans$stan)],]

# Plot it prettily
facet_labels <- found_stans %>%
  split(found_stans$feature) %>%
  sapply(function(i){paste(i$stan, i$feature, sep=": ")})
is_peak_iso %>%
  filter(feature%in%found_stans$feature) %>%
  ggplot() +
  geom_bar(aes(x=file_name, y=M_area), stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
  facet_wrap(~feature, ncol = 1, scales = "free_y",
             labeller = as_labeller(facet_labels))
ggsave(filename = paste0(pretty_folder, "internal_stan_values.pdf"), 
       device = "pdf", height = 15, width = 7.5)



# Step 1: Grab the peak areas from the pooled sample(s) & normalize to injection volume
stan_data <- lapply(found_stans$feature, function(feature_num){
  stan_name <- found_stans[found_stans$feature==feature_num, "stan"]
  is_peak_iso %>%
    filter(feature==feature_num) %>%
    mutate(stan_name=stan_name) %>%
    select(stan_name, feature, mz, rt, M_area, file_name) %>%
    left_join(bionorm_values, by="file_name") %>%
    mutate(bionorm_area=M_area/norm_vol)
})

# Step 2: compare every feature to every standard and calculate min CV
BMIS <- pbsapply(unique(complete_peaks$feature), function(feature_num){
  feature_pooled <- complete_peaks %>%
    filter(feature==feature_num) %>%
    slice(grep(pattern = "Poo", file_name)) %>%
    left_join(bionorm_values, by="file_name") %>%
    mutate(bionorm_area=M_area/norm_vol)
  initial_rsd <- sd(feature_pooled$bionorm_area)/mean(feature_pooled$bionorm_area)
  
  suppressMessages(
    stan_improvements <- stan_data %>%
      do.call(what=rbind) %>%
      slice(grep(pattern = "Poo", file_name)) %>%
      select(file_name, stan_name, stan_bionorm_area=bionorm_area) %>%
      left_join(feature_pooled, by="file_name") %>%
      mutate(stan_norm_area=bionorm_area/stan_bionorm_area) %>%
      group_by(stan_name) %>%
      summarize(MIS_rsd=sd(stan_norm_area, na.rm = TRUE)/
                  mean(stan_norm_area, na.rm=TRUE)) %>%
      ungroup() %>%
      mutate(improvement=(initial_rsd-MIS_rsd)/initial_rsd) %>%
      mutate(initial_rsd=initial_rsd) %>%
      mutate(acceptable=improvement>cut.off) %>%
      rbind(c("None", initial_rsd, 0, initial_rsd, TRUE), .) %>%
      filter(improvement==max(improvement, na.rm = TRUE))
  )
  if(nrow(stan_improvements)>1){
    return("None")
  } else {
    return(stan_improvements$stan_name)
  }
}) %>%
  data.frame(feature=names(.), BMIS=.)

# Step 3: Calculate new peak areas, normalizing to B-MIS
stan_df <- stan_data %>%
  do.call(what = rbind) %>%
  select("BMIS"=stan_name, file_name, bionorm_area) %>%
  rbind(data.frame(BMIS="None", file_name=unique(.$file_name), bionorm_area=1))
BMISed_feature_peaks <- complete_peaks %>%
  left_join(BMIS, by="feature") %>%
  left_join(stan_df, by=c("BMIS", "file_name")) %>%
  arrange(feature) %>%
  group_by(feature) %>%
  mutate(BMISed_area=(M_area/bionorm_area)*mean(bionorm_area, na.rm=TRUE)) %>%
  ungroup() %>%
  select(-bionorm_area) %>%
  filter(!feature%in%found_stans$feature)
BMISed_features <- BMISed_feature_peaks %>%
  group_by(feature) %>%
  summarise(mzmed=median(mz), rtmed=median(rt), BMIS=unique(BMIS), BMIS_avg=mean(BMISed_area, na.rm=TRUE))
write.csv(BMISed_feature_peaks, 
          file = paste0(intermediate_folder, "BMISed_feature_peaks.csv"), 
          row.names = FALSE)



# Write out peak and feature lists ----
write.csv(file = paste0(pretty_folder, "final_features.csv"),
          BMISed_features, row.names = FALSE)
write.csv(file = paste0(pretty_folder, "final_peaks.csv"), 
          BMISed_feature_peaks, row.names = FALSE)
