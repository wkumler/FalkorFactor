
# Setup things ----
library(xcms)
library(tidyverse)

load("xcms/raw_data")
x <- raw_data %>%
  filterMsLevel(msLevel. = 1L) %>%
  selectFeatureData(fcol = c(MSnbase:::.MSnExpReqFvarLabels, "centroided")) %>%
  lapply(1:length(fileNames(.)), FUN=filterFile, object = .) %>%
  `[[`(1) %>%
  spectra()


# Generate DF, pre-processing ----
mzs <- lapply(x, mz)
mz <- unlist(mzs, use.names = FALSE)
int <- unlist(lapply(x, intensity), use.names = FALSE)
rts <- unlist(lapply(x, rtime))
rt <- rep(rts, sapply(mzs, length))
all_data <- data.frame(mz, int, rt)
ppm <- 2.5


# Auxilary functions ----

peak <- setClass("peak", slots = list(roi_data="data.frame",
                                      peak_data="data.frame",
                                      start_rt="numeric",
                                      end_rt="numeric",
                                      width_rt="numeric",
                                      start_scan="numeric",
                                      end_scan="numeric",
                                      width_scan="numeric",
                                      height="numeric",
                                      area_absolute="numeric",
                                      area_above_noise="numeric",
                                      best_wavelet="numeric",
                                      ridge_length="numeric",
                                      ridge_percentage="numeric",
                                      ridge_drift="numeric",
                                      roi_sharpness="numeric",
                                      roi_accuracy="numeric",
                                      roi_noise_xcms="numeric",
                                      roi_noise_IQR="numeric",
                                      roi_noise_wopeaks="numeric",
                                      SNR_xcms="numeric",
                                      SNR_IQR="numeric",
                                      SNR_wopeaks="numeric",
                                      coef_area = "numeric"))


lmaoPlotEm <- function(roi, default_layout=T, labels = T) {
  Da_spread <- roi$mz[which.max(roi$int)]*ppm/1000000
  if(default_layout){
    layout(matrix(c(1,2), nrow = 2))
  }
  
  roi_sub_IQR <- roi$int[roi$int<median(roi$int)+IQR(roi$int)/2]
  roi_background <- median(roi_sub_IQR)
  roi_noise <- sd(roi_sub_IQR)
  
  par(mar=c(0.1, 4.1, 2.1, 0.1))
  int_colors <- hcl.colors(100, palette = "plasma")[cut(roi$int, breaks = 100)]
  plot(roi$rt, roi$mz, col=int_colors, xaxt="n", xlab="", pch=19, cex=1,
       ylim=c(min(roi$mz)*0.999999, max(roi$mz)*1.000001))
  if(labels){
    legend("topleft", legend = paste("Min m/z:", round(min(roi$mz), 5)))
    legend("topright", legend = paste("Max m/z:", round(max(roi$mz), 5)))
    legend("bottomleft", legend = paste("Actual m/z diff:", 
                                        round(max(roi$mz)-min(roi$mz), 5)))
    legend("bottomright", legend = paste("Predicted epsilon:", round(Da_spread*2, 5)))
  }
  par(mar=c(4.1, 4.1, 0.1, 0.1))
  plot(roi$rt, roi$int, col=int_colors, pch=19)
  legend("topright", legend = paste("Simple Max/Noise:", round((max(roi$int)-roi_background)/roi_noise)))
  if(default_layout){
    layout(1)
  }
}
diagnoseROI <- function(roi){
  roi_start_scan <- which(rts==roi[1, "rt"])-1
  
  # Calculate ROI "sharpness": inverse metric of signal-to-noise?
  sharpness <- 1-summary(lm(sort(roi$int)~roi$rt))$r.squared
  
  # Calculate ROI actual m/z diff vs predicted epsilon
  accuracy <- (max(roi$mz)-min(roi$mz))/
    ((roi$mz[which.max(roi$int)]*ppm/1000000)*2)
  
  
  # Wavelet transform
  # scales <- seq(1, 2^ceiling(log2(length(roi$int)))/12, length.out = 11)
  scales <- 1:(peakwidth[2]/2)
  wcoef_matrix <- xcms:::MSW.cwt(roi$int, scales, wavelet = "mexh")
  local_maxima <- xcms:::MSW.getLocalMaximumCWT(wcoef_matrix)
  possible_peaks <- xcms:::MSW.getRidge(local_maxima)
  num_scales <- length(attr(possible_peaks, "scales"))
  # Remove all peaks that have maxima in less than half the scales 
  # ADJUST THIS LATER IF NECESSARY
  possible_peaks <- possible_peaks[sapply(possible_peaks, length)>
                                     (num_scales/2)]
  # Collect ridge data on simple peaks
  ridge_lengths <- sapply(possible_peaks, length)
  ridge_percentages <- round(ridge_lengths/num_scales, digits = 2)
  ridge_drift <- sapply(possible_peaks, function(x)length(unique(x)))/ridge_lengths
  
  # Find centers by finding the scan with the highest values to left and right
  peak_center_scans <- sapply(possible_peaks, function(ridge_maxes){
    wavelet_ints <- sapply(unique(ridge_maxes), function(roi_row){
      sum(roi$int[max(1, roi_row-5):(roi_row+5)], na.rm = T)
    })
    unique(ridge_maxes)[which.max(wavelet_ints)]
  }, USE.NAMES = F)
  best_scales <- sapply(peak_center_scans, function(x){
    scale_maxes <- as.logical(local_maxima[x,])
    max(as.numeric(colnames(local_maxima))[scale_maxes])
  })
  peak_edges <- lapply(seq_along(peak_center_scans), function(x){
    left_shoulder_offset <- peak_center_scans[x]-best_scales[x]
    right_shoulder_offset <- peak_center_scans[x]+best_scales[x]
    peak_edges <- xcms:::descendMinTol(roi$int, maxDescOutlier = 2,
                         startpos = c(left_shoulder_offset, right_shoulder_offset))
  })
  peak_lefts <- sapply(peak_edges, `[`, 1)
  peak_lefts[peak_lefts<1] <- 1
  peak_rights <- sapply(peak_edges, `[`, 2)
  peak_rights[peak_rights>nrow(roi)] <- nrow(roi)
  roi_sub_IQR <- roi$int[roi$int<median(roi$int)+IQR(roi$int)/2]
  roi_noise_IQR <- c(median(roi_sub_IQR), sd(roi_sub_IQR))
  xcms_noise_baseline <- xcms:::estimateChromNoise(roi$int, trim = 0.05, 
                                                   minPts = 3*peakwidth[1])
  roi_noise_xcms <- xcms:::getLocalNoiseEstimate(roi$int, 1:nrow(roi), 
                                                 6:(nrow(roi)-5),
                                                 peakwidth*3/2, length(rts), 
                                                 xcms_noise_baseline, 8)
  roi_peak_scans <- unlist(lapply(seq_along(peak_lefts), function(x){
    seq(peak_lefts[x], peak_rights[x])}))
  roi_nonpeak_scans <- (1:nrow(roi))[-roi_peak_scans]
  if(!length(roi_nonpeak_scans)){ # If the peak runs the whole length of the ROI
    roi_nonpeak_scans <- c(1, nrow(roi)) # Use just first and last scan
  }
  roi_noise_wopeaks <- c(mean(roi$int[roi_nonpeak_scans]), 
                         sd(roi$int[roi_nonpeak_scans]))
  roi_nonpeak_scans <- (1:nrow(roi))[-roi_peak_scans]
  if(length(roi_nonpeak_scans)<2){ roi_nonpeak_scans <- c(1, nrow(roi))}
  roi_noise_wopeaks <- c(median(roi$int[roi_nonpeak_scans]), 
                         sd(roi$int[roi_nonpeak_scans]))
  
  
  layout(matrix(c(rep(1, 30), rep(2, 30), 0, rep(3, 28), 0), nrow = 3, byrow = T))
  par(mar=c(0.1, 4.1, 2.1, 0.1))
  plot(roi$rt, roi$int, type="l", lwd=1, xaxt="n", ylab="EIC intensity")
  abline(v=rts[peak_center_scans+roi_start_scan], col="red")
  abline(v=rts[peak_rights+roi_start_scan], col="blue")
  abline(v=rts[peak_lefts+roi_start_scan], col="blue")
  legend("topright", legend = c(paste("IQR SNR:", 
                                    round((max(roi$int)-roi_noise_IQR[1])/roi_noise_IQR[2])),
                                paste("xcms SNR:",
                                    round((max(roi$int)-roi_noise_xcms[1])/roi_noise_xcms[2])),
                                paste("nonpeak SNR:",
                                    round((max(roi$int)-roi_noise_wopeaks[1])/
                                            roi_noise_wopeaks[2]))))
  par(mar=c(4.1, 4.1, 0.1, 0.1))
  plot(roi$rt, wcoef_matrix[,ncol(wcoef_matrix)], 
       xlab="Retention time (s)", ylab="Wavelet coefficient", type="n")
  for(i in ncol(wcoef_matrix):1){
    lines(roi$rt, wcoef_matrix[,i], lwd=2,
          col=rainbow(ncol(wcoef_matrix))[i])
  }
  abline(h=0)
  par(mar=c(0.1, 4.1, 2.1, 0.1))
  image(wcoef_matrix, axes=F, col=hcl.colors(100, "Geyser"))
  min_scale <- as.numeric(min(colnames(wcoef_matrix)))
  yax <- pretty(colnames(wcoef_matrix))-min_scale
  axis(side = 2, at = yax/max(yax), labels = yax+min_scale, las=1, line = 1.7)
  mtext(text = "Wavelet scale", side = 2, line = 4)
  image(local_maxima, add=T, col=c("#FFFFFF00", "#000000FF"))
  layout(1)
  par(mar=c(4.1, 4.1, 0.1, 0.1))
}



# Generate ROI list ----
peakwidth <- c(20, 80)
prefilter <- c(1, 1)

data <- all_data %>% filter(mz>100&mz<120) %>% filter(rt>60&rt<1100)

pb <- txtProgressBar(min = 0, max = 1, style = 3)
roi_list <- list()
data_start_length <- nrow(data)
while(nrow(data)>0){
  setTxtProgressBar(pb, 1-(sqrt(nrow(data))/sqrt(data_start_length)))
  point_of_interest <- data[which.max(data$int),]
  epsilon_Da <- point_of_interest$mz*ppm/1000000
  upper_eic_mz <- point_of_interest$mz+epsilon_Da
  lower_eic_mz <- point_of_interest$mz-epsilon_Da
  eic <- filter(data, mz>lower_eic_mz & mz<upper_eic_mz)

  # If the ROI can't contain a peak bc too short, remove it
  if(nrow(eic)<peakwidth[1]){
    data <- data[data$mz<lower_eic_mz | data$mz>upper_eic_mz,]
    next
  }
  
  # If there aren't enough points above the prefilter threshold
  runs_above_prefilter <- rle(eic$int>prefilter[2])
  max_run_length <- max(runs_above_prefilter$lengths[runs_above_prefilter$values])
  if(max_run_length<prefilter[1]){
    # Remove the ROI
    data <- data[data$mz<lower_eic_mz | data$mz>upper_eic_mz,]
    next
  }
  
  # Split by missed scans
  roi_encoding <- rle(rts%in%eic$rt)
  peak_lengths <- roi_encoding$lengths[roi_encoding$values]
  roi_sub_list <- split(eic, rep(1:length(peak_lengths), times = peak_lengths))
  eic_split_list <- roi_sub_list[sapply(roi_sub_list, nrow)>peakwidth[1]]
  if(!length(eic_split_list)){ # If there were no reasonable ROIs found
    data <- data[data$mz<lower_eic_mz | data$mz>upper_eic_mz,]
    next
  }

  # Put it all together and save
  roi_slot <- (length(roi_list)+1):(length(roi_list)+length(eic_split_list))
  roi_list[roi_slot] <- eic_split_list
  data <- data[data$mz<lower_eic_mz | data$mz>upper_eic_mz,]
}
close(pb)





# Peak picking ----
for(i in 1:length(roi_list)){
  roi <- roi_list[[i]]
  roi_start_scan <- which(rts==roi[1, "rt"])-1
  
  # Calculate ROI "sharpness": inverse metric of signal-to-noise?
  sharpness <- 1-summary(lm(sort(roi$int)~roi$rt))$r.squared
  
  # Calculate ROI actual m/z diff vs predicted epsilon
  accuracy <- (max(roi$mz)-min(roi$mz))/
    ((roi$mz[which.max(roi$int)]*ppm/1000000)*2)
  
  
  # Wavelet transform
  # scales <- seq(1, 2^ceiling(log2(length(roi$int)))/12, length.out = 11)
  scales <- 1:(peakwidth[2]/2)
  wcoef_matrix <- xcms:::MSW.cwt(roi$int, scales, wavelet = "mexh")
  local_maxima <- xcms:::MSW.getLocalMaximumCWT(wcoef_matrix)
  possible_peaks <- xcms:::MSW.getRidge(local_maxima)
  num_scales <- length(attr(possible_peaks, "scales"))
  # Remove all peaks that have maxima in less than half the scales 
  # ADJUST THIS LATER IF NECESSARY
  possible_peaks <- possible_peaks[sapply(possible_peaks, length)>(num_scales/2)]
  
  # Collect ridge data on peaks
  ridge_lengths <- sapply(possible_peaks, length)
  ridge_percentages <- round(ridge_lengths/num_scales, digits = 2)
  ridge_drift <- sapply(possible_peaks, function(x)length(unique(x)))/ridge_lengths
  
  # Find centers by finding the scan with the highest values to left and right
  peak_center_scans <- sapply(possible_peaks, function(ridge_maxes){
    wavelet_ints <- sapply(unique(ridge_maxes), function(roi_row){
      sum(roi$int[max(1, roi_row-5):(roi_row+5)], na.rm = T)
    })
    unique(ridge_maxes)[which.max(wavelet_ints)]
  }, USE.NAMES = F)
  
  # Find the best wavelet scale for each peak 
  # I.e. the largest one whose maxima matched data max
  best_scales <- sapply(peak_center_scans, function(x){
    scale_maxes <- as.logical(local_maxima[x,])
    max(as.numeric(colnames(local_maxima))[scale_maxes])
  })
  
  # Calculate peak edges xcms way
  peak_edges <- lapply(seq_along(peak_center_scans), function(x){
    left_shoulder_offset <- peak_center_scans[x]-best_scales[x]
    right_shoulder_offset <- peak_center_scans[x]+best_scales[x]
    xcms:::descendMinTol(roi$int, maxDescOutlier = 2,
                         startpos = c(left_shoulder_offset, right_shoulder_offset))
  })
  peak_lefts <- sapply(peak_edges, `[`, 1)
  peak_lefts[peak_lefts<1] <- 1 # Make sure they're all on the map
  peak_rights <- sapply(peak_edges, `[`, 2)
  peak_rights[peak_rights>nrow(roi)] <- nrow(roi) # Same as above
  
  # Calculate noise for the ROI as a whole
  # (IQR method)
  roi_sub_IQR <- roi$int[roi$int<median(roi$int)+IQR(roi$int)]
  roi_noise_IQR <- c(median(roi_sub_IQR), sd(roi_sub_IQR))
  # (xcms method)
  xcms_noise_baseline <- xcms:::estimateChromNoise(roi$int, trim = 0.05, 
                                                   minPts = 3*peakwidth[1])
  roi_noise_xcms <- xcms:::getLocalNoiseEstimate(roi$int, 1:nrow(roi), 
                                                 6:(nrow(roi)-5),
                                                 peakwidth*3/2, length(rts), 
                                                 xcms_noise_baseline, 8)
  if(any(roi_noise_xcms==1)){warning("Peak has a noise value of 1")}
  # (peak removal method)
  roi_peak_scans <- unlist(lapply(seq_along(peak_lefts), function(x){
    seq(peak_lefts[x], peak_rights[x])}))
  roi_nonpeak_scans <- (1:nrow(roi))[-roi_peak_scans]
  if(length(roi_nonpeak_scans)<2){ # If the peak runs the whole length of the ROI
    roi_nonpeak_scans <- c(1, nrow(roi)) # Use just first and last scan
  }
  roi_noise_wopeaks <- c(median(roi$int[roi_nonpeak_scans]), 
                         sd(roi$int[roi_nonpeak_scans]))
  
  # Calculate absolute and relative peak area
  peak_areas <- lapply(seq_along(peak_center_scans), function(x){
    given_peak_scans <- peak_lefts[x]:peak_rights[x]
    idx = 2:length(given_peak_scans)
    y <- roi$int[given_peak_scans]
    absolute_area <- ((given_peak_scans[idx] - given_peak_scans[idx-1]) %*% 
                        (y[idx] + y[idx-1]))/2
    
    area_above_noise <- absolute_area-roi_noise_xcms[x]*length(given_peak_scans)
    
    return(c(absolute_area, area_above_noise))
  })
  absolute_areas <- sapply(peak_areas, `[[`, 1)
  relative_areas <- sapply(peak_areas, `[[`, 2)
  
  # Calculate best coefficient:area estimate
  coef_areas <- sapply(seq_along(peak_areas), function(x){
    given_peak_scans <- peak_lefts[x]:peak_rights[x]
    given_peak_coefs <- wcoef_matrix[given_peak_scans, ]
    return(max(given_peak_coefs)/absolute_areas[x])
  })
}



# Post-processing ----

lmaoPlotEm(roi_list[[1]])
diagnoseWavelets(roi_list[[1]])

roi <- roi_list[[sample(length(roi_list), 1)]]
diagnoseROI(roi)
lmaoPlotEm(roi)

