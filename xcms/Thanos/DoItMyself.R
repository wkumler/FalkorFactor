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

ROI <- setClass("ROI", slots = list(data="data.frame", 
                                    ROI_number="numeric",
                                    mz_info="list",
                                    consistency="rle"))



# Generating DF, pre-processing ----
mzs <- lapply(x, mz)
mz <- unlist(mzs, use.names = FALSE)
int <- unlist(lapply(x, intensity), use.names = FALSE)
rts <- unlist(lapply(x, rtime))
rt <- rep(rts, sapply(mzs, length))
all_data <- data.frame(mz, int, rt)
ppm <- 2.5

# Optionally, remove all orphan data points
#seq_data <- subset(all_data, duplicated(mz) | duplicated(mz, fromLast=TRUE))

# Determine typical m/z spread for a couple good peaks
# Unnecessary because m/z will be determined on-the-fly by ROI calculation, but looks nice
lmaoPlotEm <- function(data_i, default_layout=T, labels = T) {
  Da_spread <- data_i$mz[which.max(data_i$int)]*ppm/1000000
  if(default_layout){
    layout(matrix(c(1,2), nrow = 2))
  }
  par(mar=c(0.1, 4.1, 2.1, 0.1))
  int_colors <- hcl.colors(100, palette = "plasma")[cut(data_i$int, breaks = 100)]
  plot(data_i$rt, data_i$mz, col=int_colors, xaxt="n", xlab="", pch=19, cex=1,
       ylim=c(min(data_i$mz)*0.999999, max(data_i$mz)*1.000001))
  if(labels){
    legend("topleft", legend = paste("Min m/z:", round(min(data_i$mz), 5)))
    legend("topright", legend = paste("Max m/z:", round(max(data_i$mz), 5)))
    legend("bottomleft", legend = paste("Actual m/z diff:", 
                                        round(max(data_i$mz)-min(data_i$mz), 5)))
    legend("bottomright", legend = paste("Predicted epsilon:", round(Da_spread*2, 5)))
  }
  par(mar=c(4.1, 4.1, 0.1, 0.1))
  plot(data_i$rt, data_i$int, col=int_colors, pch=19)
  if(default_layout){
    layout(1)
  }
}
# all_data %>% filter(mz>90.0910&mz<90.0925) %>% lmaoPlotEm()
# all_data %>% filter(mz>100.023&mz<100.025) %>% lmaoPlotEm()
# all_data %>% filter(mz>117.048&mz<117.049) %>% lmaoPlotEm()
# all_data %>% filter(mz>285.040&mz<285.050) %>% lmaoPlotEm()
# all_data %>% filter(mz>800.810&mz<800.815) %>% lmaoPlotEm()
# all_data %>% filter(mz>135.0&mz<135.1) %>% lmaoPlotEm()
# all_data %>% filter(mz>96.0&mz<96.1) %>% lmaoPlotEm()



# Generate ROI list ----
# Sort by intensity
ppm <- 2.5
peakwidth <- c(20, 80)
prefilter <- c(3, 100)

roi_list <- list()

data <- all_data %>% filter(mz>100&mz<120)

pb <- txtProgressBar(min = 0, max = 1, style = 3)
roi_start_length <- nrow(data)
while(nrow(data)>0){
  setTxtProgressBar(pb, 1-(sqrt(nrow(data))/sqrt(roi_start_length)))
  point_of_interest <- data[which.max(data$int),]
  epsilon_Da <- point_of_interest$mz*ppm/1000000
  upper_roi_mz <- point_of_interest$mz+epsilon_Da
  lower_roi_mz <- point_of_interest$mz-epsilon_Da
  roi <- filter(data, mz>lower_roi_mz & mz<upper_roi_mz)
  #QUESTION: Should the ROIs be sliced when they miss scans? Or drop below an intensity threshold?
  
  # If the ROI can't contain a peak bc too short, remove it
  if(nrow(roi)<peakwidth[1]){
    data <- data[data$mz<lower_roi_mz | data$mz>upper_roi_mz,]
    next
  }
  
  # If there aren't enough points above the prefilter threshold
  runs_above_prefilter <- rle(roi$int>prefilter[2])
  max_run_length <- max(runs_above_prefilter$lengths[runs_above_prefilter$values])
  if(max_run_length<prefilter[1]){
    # Remove the ROI
    data <- data[data$mz<lower_roi_mz | data$mz>upper_roi_mz,]
    next
  }
  
  # Calculate ROI consistency by checking for missed scans
  consistency <- rle(rts%in%roi$rt)
  
  # Calculate ROI "sharpness": inverse metric of signal-to-noise?
  sharpness <- 1/summary(lm(sort(roi$int)~roi$rt))$r.squared
  
  # Other ROI-wide information can go here
  
  
  # Wavelet transform
  # scales <- seq(peakwidth[1]/2, peakwidth[2]/2)
  # Scales now run from 1+ to catch sharp peaks w/o surrounding noise
  # See peak at mz>119.99&mz<119.994 which is missed because EIC is too short to support wavelets 10+
  scales <- seq(1, 2^ceiling(log2(length(roi$int)))/12, length.out = 11)
  scales <- 1:20
  # QUESTION: since only intensity is passed to the wavelet function, how does it handle missed scans?
  possible_peaks <- xcms:::MSW.cwt(roi$int, scales, wavelet = "mexh") %>%
    xcms:::MSW.getLocalMaximumCWT() %>%
    xcms:::MSW.getRidge()
  peak_centers <- sapply(possible_peaks, function(ridge_maxes){
    wavelet_ints <- sapply(unique(ridge_maxes), function(roi_row){
      sum(roi$int[max(1, roi_row-5):(roi_row+5)], na.rm = T)
    })
    unique(ridge_maxes)[which.max(wavelet_ints)]
  })
  ridge_lengths <- sapply(possible_peaks, length)
  ridge_percentages <- round(ridge_lengths/length(attr(possible_peaks, "scales")),
                             digits = 2)
  ridge_drift <- sapply(possible_peaks, function(x)length(unique(x)))/ridge_lengths
  
  # Put it all together and save
  roi_list[[length(roi_list)+1]] <- 
    ROI(data = roi,
        ROI_number = length(roi_list)+1,
        mz_info = list("wmean"=weighted.mean(roi$mz, roi$int),
                       "mz_min"=min(roi$mz), "mz_max"=max(roi$mz)),
        consistency = consistency)
  
  data <- data[data$mz<lower_roi_mz | data$mz>upper_roi_mz,]
}
close(pb)


# Post-processing ----


diagnoseWavelets <- function(roi){
  wcoef_matrix <- xcms:::MSW.cwt(roi$int, scales = scales, wavelet = "mexh")
  local_maxima <- xcms:::MSW.getLocalMaximumCWT(wcoef_matrix)
  
  layout(matrix(c(rep(1, 30), rep(2, 30), 0, rep(3, 28), 0), nrow = 3, byrow = T))
  par(mar=c(0.1, 4.1, 2.1, 0.1))
  plot(roi$rt, roi$int, type="l", lwd=1, xaxt="n", ylab="EIC intensity")
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
}
