# williWave_source
# Will re-writing the CentWave algorithm


# Setup things ----
library(xcms)
library(dplyr)

load("xcms/raw_data")
x <- raw_data %>%
  filterMsLevel(msLevel. = 1L) %>%
  selectFeatureData(fcol = c(MSnbase:::.MSnExpReqFvarLabels, "centroided")) %>%
  lapply(1:length(fileNames(.)), FUN=filterFile, object = .) %>%
  `[[`(1) %>%
  spectra()


# Generate ROIs ----
mzs <- lapply(x, mz)
val_count <- cumsum(lengths(mzs, FALSE))
mz <- unlist(mzs, use.names = FALSE)
int <- unlist(lapply(x, intensity), use.names = FALSE)
rt <- unlist(lapply(x, rtime), use.names = FALSE)
rts <- rep(rt, sapply(mzs, length))

scanindex <- as.integer(c(0, val_count[-length(val_count)]))
scantime = rt
scanrange <- c(1, length(scantime))
mz_span <- c(0.0005) # Maximum spread of m/z values across a well-defined peak, plus some buffer
ppm <- ceiling((mz_span*1000000)/132)
peakwidth <- c(20, 80)
min_peak_width <- min(peakwidth)/2
#min_centroids <- max(4, min_peak_width - 2)
min_centroids <- 9 #To parallel milliWave
prefilter = c(3, 10000)
noise = 0

roi_list <- .Call("findmzROI", 
                 mz, int, scanindex, 
                 as.double(c(0, 0)), 
                 as.integer(scanrange), 
                 as.integer(length(scantime)), 
                 as.double(ppm * 1e-06), 
                 as.integer(min_centroids), 
                 as.integer(prefilter), 
                 as.integer(noise), PACKAGE = "xcms")


# For each ROI ----
# Get EIC
f <- 500
roi <- roi_list[[f]]
mzrange <- c(roi$mzmin, roi$mzmax)
eic_span <- c(max(scanrange[1], roi$scmin - max(peakwidth)*3/2), 
              min(scanrange[2], roi$scmax + max(peakwidth)*3/2))
eic <- .Call("getEIC", mz, int, scanindex, as.double(mzrange), 
             as.integer(eic_span), as.integer(length(scanindex)), PACKAGE = "xcms")
eic_scan_start <- min(eic$scan)


# Get wavelets
#scales <- seq.int(peakwidth[1]/2, peakwidth[2]/2)
scales <- 11:44 #To parallel milliWave
w_coefs <- xcms:::MSW.cwt(eic$intensity, scales = scales, wavelet = "mexh")


# Get Ridgelines
local_maxima <- xcms:::MSW.getLocalMaximumCWT(w_coefs)
ridgelines <- xcms:::MSW.getRidge(local_maxima)
ridgeline_maxima <- lapply(ridgelines, function(x){x+eic_scan_start})
ridge_length <- sapply(ridgelines, length)
ridge_percentage <- round(ridge_length/ncol(w_coefs), digits = 2)


# Find all local maxima in wavelets
peaks_in_roi <- matrix(numeric(6*length(ridgeline_maxima)), ncol = 6, 
                       dimnames = list(NULL, c("best_wavelet_scale", 
                                               "peak_midpoint", "peak_left_edge",
                                               "peak_right_edge", "peak_height",
                                               "ridge_intensity")))
for(wavelet_peak in seq_along(ridgeline_maxima)){
  # Calculate peak area around the wavelet peak to find the best fitting wavelet for each peak
  wavelet_ints <- sapply(unique(ridgeline_maxima[[wavelet_peak]]), function(x){
    left_bound <- max(eic_scan_start, x-min_peak_width/2)
    right_bound <- min(max(eic$scan), x+min_peak_width/2)
    sum(eic$intensity[(left_bound:right_bound)-eic_scan_start])
  })
  
  # Find the wavelet scans with the best match
  best_wavelet_scan <- unique(ridgeline_maxima[[wavelet_peak]])[which.max(wavelet_ints)]
  # Find the wavelet scales with the best match
  best_wavelets <- scales[as.logical(local_maxima[as.numeric(best_wavelet_scan)-eic_scan_start,])]
  best_wavelet <- max(best_wavelets)
  
  # Find peak edges
  left_shoulder_offset <- (best_wavelet_scan-best_wavelet)-eic_scan_start
  right_shoulder_offset <- (best_wavelet_scan+best_wavelet)-eic_scan_start
  peak_edges <- xcms:::descendMinTol(eic$intensity, 
                                     startpos = c(left_shoulder_offset, right_shoulder_offset),
                                     maxDescOutlier = min_peak_width)+eic_scan_start
  # Make sure the edges make sense
  if(peak_edges[1]<eic$scan[1]){
    peak_edges[1] <- eic$scan[1]
  }
  if(peak_edges[2]>eic$scan[length(eic$scan)]){
    peak_edges[2] <- eic$scan[length(eic$scan)]
  }
  
  # Find the height of the peak
  peak_height <- max(eic$intensity[(peak_edges[1]:peak_edges[2])-eic_scan_start])
  
  # Estimate peak noise (xcms algorithm)
  xcms_noise_baseline <- xcms:::estimateChromNoise(eic$intensity, trim = 0.05, minPts = 3*min_peak_width)
  peak_plus_noise <- max(eic$scan[1], roi$scmin - floor(min(peakwidth)/2)):
    min(eic$scan[length(eic$scan)], roi$scmax + floor(min(peakwidth)/2))
  noiserange <- c(min(scales) * 3, max(scales) * 3)
  minPtsAboveBaseLine <- max(4, min_peak_width - 2)
  lnoise <- xcms:::getLocalNoiseEstimate(eic$intensity, eic$scan, peak_plus_noise, 
                                         noiserange, length(scantime), 
                                         threshold = xcms_noise_baseline, 
                                         num = minPtsAboveBaseLine)
  # Estimate peak noise (ADAP algorithm)
  ADAP_noise <- estimateADAPnoise(eic$scan, eic$intensity, peak_edges)
  
  # 
  
  # Collect all the information
  peakinfo <- c("best_wavelet_scale"=best_wavelet,
                "peak_midpoint"=best_wavelet_scan,
                "peak_left_edge"=peak_edges[1],
                "peak_right_edge"=peak_edges[2],
                "peak_height"=peak_height,
                "ridge_intensity"=sum(w_coefs[best_wavelet_scan, ])/ncol(w_coefs))
  peaks_in_roi[wavelet_peak,] <- peakinfo
}
peak_table_roi <- cbind(peaks_in_roi, ridge_length, ridge_percentage)




# Calculate peak area #NEEDS BASELINE ESTIMATE
peak_idx <- (peak_edges[1]-eic_scan_start):(peak_edges[2]-eic_scan_start)
peak_scans <- eic$scan[peak_idx]
peak_ints <- eic$intensity[peak_idx]
riemann_vals <- 2:length(peak_scans)
integral <- as.double((peak_scans[riemann_vals] - peak_scans[riemann_vals-1]) %*% 
                        (peak_ints[riemann_vals] + peak_ints[riemann_vals-1]))/2
# Debug peak area calculation
# xleft = peak_edges[1]
# xright = peak_edges[2]
# height = integral/(xright-xleft)
# lines(eic$scan, eic$intensity, lwd=2)
# rect(xleft, 0, xright, height, lwd=2)


# Filters ----
# getPointsAboveThreshold()
# getSignal2Noise()
# getWaveletCoef1()
# getWaveletCoef2()
# getRidgelineIntensity()
# getRidgelineLength()
# getCoef2Area()
# getPeakShape()

# Diagnostics ----
diagnoseROI <- function(ROI_number){
  roi <- roi_list[[ROI_number]]
  mzrange <- c(roi$mzmin, roi$mzmax)
  eic_span <- c(max(scanrange[1], roi$scmin - max(peakwidth)*3/2), 
                min(scanrange[2], roi$scmax + max(peakwidth)*3/2))
  eic <- .Call("getEIC", mz, int, scanindex, as.double(mzrange), 
               as.integer(eic_span), as.integer(length(scanindex)), PACKAGE = "xcms")
  eic_scan_start <- min(eic$scan)
  
  
  # Get wavelets
  #scales <- seq.int(peakwidth[1]/2, peakwidth[2]/2)
  scales <- 11:44 #To parallel milliWave
  w_coefs <- xcms:::MSW.cwt(eic$intensity, scales = scales, wavelet = "mexh")
  
  
  
  # Get Ridgelines
  local_maxima <- xcms:::MSW.getLocalMaximumCWT(w_coefs)
  ridgelines <- xcms:::MSW.getRidge(local_maxima)
  ridgeline_maxima <- lapply(ridgelines, function(x){x+eic_scan_start})
  ridge_length <- sapply(ridgelines, length)
  ridge_percentage <- ridge_length/length(scales)
  
  # Find all local maxima in wavelets
  # Preallocate matrix
  peaks_in_roi <- matrix(numeric(5*length(ridgeline_maxima)), ncol = 5, 
                         dimnames = list(NULL, c("best_wavelet_scale", 
                                                 "peak_midpoint", "peak_left_edge",
                                                 "peak_right_edge", "peak_height")))
  for(wavelet_peak in seq_along(ridgeline_maxima)){
    wavelet_ints <- sapply(unique(ridgeline_maxima[[wavelet_peak]]), function(x){
      # Calculate peak area around the wavelet peak to find the best fitting wavelet
      left_bound <- max(eic_scan_start, x-min_peak_width/2)
      right_bound <- min(max(eic$scan), x+min_peak_width/2)
      sum(eic$intensity[(left_bound:right_bound)-eic_scan_start])
    })
    # Find the wavelet scans with the best match
    best_wavelet_scan <- unique(ridgeline_maxima[[wavelet_peak]])[which.max(wavelet_ints)]
    # Find the wavelet scales with the best match
    best_wavelets <- scales[as.logical(local_maxima[as.numeric(best_wavelet_scan)-eic_scan_start,])]
    best_wavelet <- max(best_wavelets)
    # Find peak edges
    left_shoulder_offset <- (best_wavelet_scan-best_wavelet)-eic_scan_start
    right_shoulder_offset <- (best_wavelet_scan+best_wavelet)-eic_scan_start
    peak_edges <- xcms:::descendMinTol(eic$intensity, 
                                       startpos = c(left_shoulder_offset, right_shoulder_offset),
                                       maxDescOutlier = min_peak_width)+eic_scan_start
    if(peak_edges[1]<eic$scan[1]){
      # warning(paste0("Found a peak in ROI ",f, ", wavelet peak ", 
      #                wavelet_peak, ", off the left edge of the map!"))
      peak_edges[1] <- eic$scan[1]
    }
    if(peak_edges[2]>eic$scan[length(eic$scan)]){
      # warning(paste0("Found a peak in ROI ",f, ", wavelet peak ", 
      #                wavelet_peak, ", off the right edge of the map!"))
      peak_edges[2] <- eic$scan[length(eic$scan)]
    }
    peak_height <- max(eic$intensity[(peak_edges[1]:peak_edges[2])-eic_scan_start])
    peakinfo <- c("best_wavelet_scale"=best_wavelet,
                  "peak_midpoint"=best_wavelet_scan,
                  "peak_left_edge"=peak_edges[1],
                  "peak_right_edge"=peak_edges[2],
                  "peak_height"=peak_height)
    peaks_in_roi[wavelet_peak,] <- peakinfo
  }
  peak_table_roi <- cbind(peaks_in_roi, ridge_length, ridge_percentage)
  
  layout(matrix(c(rep(1, 30), rep(2, 30), 0, rep(3, 28), 0), nrow = 3, byrow = T))
  par(mar=c(0.1, 4.1, 2.1, 0.1))
  plot(eic$scan, eic$intensity, type="l", lwd=1, xaxt="n", ylab="EIC intensity",
       main=paste0("Peak at ", roi$mz, " m/z, ROI: ", ROI_number))
  par(mar=c(0.1, 4.1, 0.1, 0.1))
  plot(eic$scan, w_coefs[,dim(w_coefs)[2]], xaxt="n", ylab="Wavelet coefficient")
  for(i in dim(w_coefs)[2]:1){
    points(eic$scan, w_coefs[,i], col=rainbow(dim(w_coefs)[2])[i], cex=0.5, pch=19)
  }
  abline(h=0, lwd=2)
  par(mar=c(4.1, 4.1, 0.1, 0.1))
  image(w_coefs, axes=F, col=hcl.colors(100, "Geyser"))
  xax <- pretty(eic$scan)-min(pretty(eic$scan))
  min_scale <- as.numeric(min(colnames(w_coefs)))
  yax <- pretty(colnames(w_coefs))-min_scale
  axis(side = 1, at = xax/max(xax), labels = xax+min(pretty(eic$scan)))
  axis(side = 2, at = yax/max(yax), labels = yax+min_scale, las=1, line = 1.7)
  mtext(text = "Retention time (s)", side = 1, line = 3)
  mtext(text = "Wavelet scale", side = 2, line = 4)
  image(local_maxima, add=T, col=c("#FFFFFF00", "#000000FF"))
  layout(1)
  
  return(peak_table_roi)
}


# Additional functions ----
estimateADAPnoise <- function(eic_scans, eic_intensity, peak_edges){
  peak_width <- peak_edges[2]-peak_edges[1]
  left_bound <- max(min(eic$scan), peak_edges[1]-8*peak_width)
  right_bound <- min(max(eic$scan), peak_edges[2]+8*peak_width)
  
  window_list <- list()
  window_offset <- 0
  window <- numeric(peak_width*2)
  while(length(window)>peak_width){
    window <- c(((peak_edges[1]-peak_width):peak_edges[1])-window_offset, 
                window_offset+(peak_edges[2]:(peak_edges[2]+peak_width)))
    window <- window[window>=left_bound&window<=right_bound]
    window_offset <- window_offset+1
    window_list[[window_offset]] <- eic$intensity[window]
  }
  method1_min_sd <- min(sapply(window_list, sd))
  
  right_window_idxs <- (peak_edges[2]+peak_width):right_bound
  right_window_idxs <- right_window_idxs[right_window_idxs<right_bound]
  windows_right <- lapply(right_window_idxs, function(x){eic$intensity[peak_edges[2]:x]})
  left_window_idxs <- left_bound:(peak_edges[1]-peak_width)
  left_window_idxs <- left_window_idxs[left_window_idxs>0]
  windows_left <- lapply(left_window_idxs, function(x){eic$intensity[x:peak_edges[2]]})
  if(length(windows_left)<length(windows_right)){
    windows_left <- c(windows_left, rep(windows_left[length(windows_left)], 
                                        length(windows_right)-length(windows_left)))
  } else {
    windows_right <- c(windows_right, rep(windows_right[length(windows_right)], 
                                        length(windows_left)-length(windows_right)))
  }
  windows_increasing <- lapply(seq_along(windows_right), 
                               function(x){c(windows_right[[x]], 
                                             windows_left[[x]])})
  windows_increasing_sd <- min(sapply(windows_increasing, sd))
  
  windows_right <- lapply(right_window_idxs, function(x){eic$intensity[x:right_bound]})
  windows_left <- lapply(left_window_idxs, function(x){eic$intensity[left_bound:x]})
  if(length(windows_left)<length(windows_right)){
    windows_left <- c(windows_left, rep(windows_left[length(windows_left)], 
                                        length(windows_right)-length(windows_left)))
  } else {
    windows_right <- c(windows_right, rep(windows_right[length(windows_right)], 
                                          length(windows_left)-length(windows_right)))
  }
  windows_decreasing <- lapply(seq_along(windows_right), 
                               function(x){c(windows_right[[x]], 
                                             windows_left[[x]])})
  windows_decreasing_sd <- min(sapply(windows_decreasing, sd))
  return(min(method1_min_sd, windows_increasing_sd, windows_decreasing_sd))
}

estimateIQRnoise <- function(eic_intensity){
  
}
