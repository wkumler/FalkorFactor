
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
splitByMissedScan <- function(roi){
  roi_encoding <- rle(rts%in%roi$rt)
  if(roi_encoding$values[1]){
    peak_lengths <- roi_encoding$lengths[c(T,F)]
  } else {
    peak_lengths <- roi_encoding$lengths[c(F,T)]
  }
  roi_sub_list <- split(roi, rep(1:length(peak_lengths), times = peak_lengths))
  # Only keep segments that might contain a peak
  roi_sub_list <- roi_sub_list[sapply(roi_sub_list, nrow)>peakwidth[1]]
}

lmaoPlotEm <- function(eic, default_layout=T, labels = T) {
  Da_spread <- eic$mz[which.max(eic$int)]*ppm/1000000
  if(default_layout){
    layout(matrix(c(1,2), nrow = 2))
  }
  
  roi_sub_IQR <- roi$int[roi$int<median(roi$int)+IQR(roi$int)/2]
  roi_background <- median(roi_sub_IQR)
  roi_noise <- sd(roi_sub_IQR)
  
  par(mar=c(0.1, 4.1, 2.1, 0.1))
  int_colors <- hcl.colors(100, palette = "plasma")[cut(eic$int, breaks = 100)]
  plot(eic$rt, eic$mz, col=int_colors, xaxt="n", xlab="", pch=19, cex=1,
       ylim=c(min(eic$mz)*0.999999, max(eic$mz)*1.000001))
  if(labels){
    legend("topleft", legend = paste("Min m/z:", round(min(eic$mz), 5)))
    legend("topright", legend = paste("Max m/z:", round(max(eic$mz), 5)))
    legend("bottomleft", legend = paste("Actual m/z diff:", 
                                        round(max(eic$mz)-min(eic$mz), 5)))
    legend("bottomright", legend = paste("Predicted epsilon:", round(Da_spread*2, 5)))
  }
  par(mar=c(4.1, 4.1, 0.1, 0.1))
  plot(eic$rt, eic$int, col=int_colors, pch=19)
  legend("topright", legend = paste("Simple Max/Noise:", round((max(roi$int)-roi_background)/roi_noise)))
  if(default_layout){
    layout(1)
  }
}
diagnoseWavelets <- function(roi){
  wcoef_matrix <- xcms:::MSW.cwt(roi$int, scales = scales, wavelet = "mexh")
  local_maxima <- xcms:::MSW.getLocalMaximumCWT(wcoef_matrix)
  
  layout(matrix(c(rep(1, 30), rep(2, 30), 0, rep(3, 28), 0), nrow = 3, byrow = T))
  par(mar=c(0.1, 4.1, 2.1, 0.1))
  plot(roi$rt, roi$int, type="l", lwd=1, xaxt="n", ylab="EIC intensity")
  possible_peaks <- xcms:::MSW.getRidge(local_maxima)
  possible_peaks <- possible_peaks[sapply(possible_peaks, length)>
                                     (length(attr(possible_peaks, "scales"))/2)]
  peak_centers <- rts[sapply(possible_peaks, function(ridge_maxes){
    wavelet_ints <- sapply(unique(ridge_maxes), function(roi_row){
      sum(roi_data$int[max(1, roi_row-5):(roi_row+5)], na.rm = T)
    })
    unique(ridge_maxes)[which.max(wavelet_ints)]
  })+which(rts==roi[1, "rt"])-1]
  abline(v=peak_centers)
  
  
  par(mar=c(4.1, 4.1, 0.1, 0.1))
  plot(roi$rt, wcoef_matrix[,ncol(wcoef_matrix)], 
       xlab="Retention time (s)", ylab="Wavelet coefficient", type="n")
  for(i in ncol(wcoef_matrix):1){
    lines(roi$rt, wcoef_matrix[,i], lwd=2,
          col=rainbow(ncol(wcoef_matrix))[i])
  }
  abline(h=0, lwd=2)
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
  eic_split_list <- splitByMissedScan(eic)
  if(!length(eic_split_list)){ # If there were no reasonable ROIs found
    data <- data[data$mz<lower_eic_mz | data$mz>upper_eic_mz,]
    next
  }

  # Put it all together and save
  roi_list[(length(roi_list)+1):(length(roi_list)+length(eic_split_list))] <- eic_split_list
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
  
  # Calculate ROI noise (IQR version)
  roi_sub_IQR <- roi$int[roi$int<IQR(roi$int)]
  roi_background <- median(roi_sub_IQR)
  roi_noise <- sd(roi_sub_IQR)
  
  
  # Wavelet transform
  # scales <- seq(1, 2^ceiling(log2(length(roi$int)))/12, length.out = 11)
  scales <- 1:(peakwidth[2]/2)
  wcoef_matrix <- xcms:::MSW.cwt(roi$int, scales, wavelet = "mexh")
  possible_peaks <- xcms:::MSW.getLocalMaximumCWT(wcoef_matrix) %>%
    xcms:::MSW.getRidge()
  num_scales <- length(attr(possible_peaks, "scales"))
  # Remove all peaks that have maxima in less than half the scales ADJUST THIS LATER IF NECESSARY
  possible_peaks <- possible_peaks[sapply(possible_peaks, length)>
                                     (num_scales/2)]
  # Collect ridge data on simple peaks
  ridge_lengths <- sapply(possible_peaks, length)
  ridge_percentages <- round(ridge_lengths/num_scales, digits = 2)
  ridge_drift <- sapply(possible_peaks, function(x)length(unique(x)))/ridge_lengths
  
  # Find centers by finding the scan with the highest values to left and right
  peak_centers <- rts[sapply(possible_peaks, function(ridge_maxes){
    wavelet_ints <- sapply(unique(ridge_maxes), function(roi_row){
      sum(roi$int[max(1, roi_row-5):(roi_row+5)], na.rm = T)
    })
    unique(ridge_maxes)[which.max(wavelet_ints)]
  })+roi_start_scan]
  
  # Calculate peak edges
  # left_shoulder_offset <-
  # peak_edges <- xcms:::descendMinTol(peak_data$int, maxDescOutlier = min_peak_width,
  #                                    startpos = c(left_shoulder_offset,
  #                                                 right_shoulder_offset))
  # Calculate ROI noise (IQR method)
  roi_sub_IQR <- roi$int[roi$int<median(roi$int)+IQR(roi$int)/2]
  roi_background_IQR <- median(roi_sub_IQR)
  roi_noise_IQR <- sd(roi_sub_IQR)
  # Calculate peak noise (xcms method)
  
  # Calculate peak noise (peak removal method)
}



# Post-processing ----

lmaoPlotEm(roi_list[[1]])
diagnoseWavelets(roi_list[[1]])

roi <- roi_list[[sample(length(roi_list), 1)]]
lmaoPlotEm(roi)
diagnoseWavelets(roi)