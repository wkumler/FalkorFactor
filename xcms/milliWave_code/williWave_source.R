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


scanindex <- as.integer(c(0, val_count[-length(val_count)]))
scantime = rt
scanrange <- c(1, length(scantime))
mz_span <- c(0.0005) # Maximum spread of m/z values across a well-defined peak, plus some buffer
ppm <- ceiling((mz_span*1000000)/132)
peakwidth <- c(20, 80)
min_peak_width <- min(peakwidth)/2
min_centroids <- max(4, min_peak_width - 2)
prefilter = c(3, 100000)
noise = 0

roiList <- .Call("findmzROI", 
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
roi <- roiList[[f]]
mzrange <- c(roi$mzmin, roi$mzmax)
eic_span <- c(max(scanrange[1], roi$scmin - max(peakwidth)*3/2), 
              min(scanrange[2], roi$scmax + max(peakwidth)*3/2))
eic <- .Call("getEIC", mz, int, scanindex, as.double(mzrange), 
             as.integer(eic_span), as.integer(length(scanindex)), PACKAGE = "xcms")
eic_scan_start <- min(eic$scan)
plot(eic$scan, eic$intensity, type="l", lwd=2, ylim=c(-max(eic$intensity), 2*max(eic$intensity)))


# Get wavelets
scales <- seq.int(peakwidth[1]/2, peakwidth[2]/2)
wCoefs <- xcms:::MSW.cwt(eic$intensity, scales = scales, wavelet = "mexh")
for(i in 1:dim(wCoefs)[2]){
  points(eic$scan, wCoefs[,i], col=rainbow(dim(wCoefs)[2])[i], cex=0.5, pch=19)
}


# Get Ridgelines
local_maxima <- xcms:::MSW.getLocalMaximumCWT(wCoefs)
ridgelines <- xcms:::MSW.getRidge(local_maxima)
ridgeline_by_scan <- lapply(ridgelines, function(x){x+eic_scan_start})
### Throw out all ridges that aren't even close to the peak - ASSESS THIS LATER
good_ridgeline <- sapply(ridgeline_by_scan, function(x){
  median(x)>roi$scmin&median(x)<roi$scmax
  })
ridgeline_by_scan <- ridgeline_by_scan[good_ridgeline][[1]]
wavelet_ints <- sapply(unique(ridgeline_by_scan), function(x){
  # Calculate peak area around the wavelet peak to find the best fitting wavelet
  sum(eic$intensity[((x-min_peak_width/2):(x+min_peak_width/2))-eic_scan_start])
})
# Find the wavelet scans with the best match
best_wavelet_scan <- unique(ridgeline_by_scan)[which.max(wavelet_ints)]
# Find the wavelet scales with the best match
best_wavelets <- scales[as.logical(local_maxima[best_wavelet_scan-eic_scan_start,])]
best_wavelet <- max(best_wavelets)


# Find peak edges
left_shoulder_offset <- best_wavelet_scan-best_wavelet-eic_scan_start
right_shoulder_offset <- best_wavelet_scan+best_wavelet-eic_scan_start
peak_edges <- xcms:::descendMinTol(eic$intensity, 
                     startpos = c(left_shoulder_offset, right_shoulder_offset),
                     maxDescOutlier = min_peak_width)+eic_scan_start


# Calculate peak area
peak_idx <- (peak_edges[1]-eic_scan_start):(peak_edges[2]-eic_scan_start)
peak_scans <- eic$scan[peak_idx]
peak_ints <- eic$intensity[peak_idx]
riemann_vals <- 2:length(peak_scans)
integral <- as.double((peak_scans[riemann_vals] - peak_scans[riemann_vals-1]) %*% 
                         (peak_ints[riemann_vals] + peak_ints[riemann_vals-1]))/2

xleft = peak_edges[1]
xright = peak_edges[2]
height = integral/(xright-xleft)
lines(eic$scan, eic$intensity, lwd=2)
rect(xleft, 0, xright, height, lwd=2)


# Filters ----
getPointsAboveThreshold()
getSignal2Noise()
getWaveletCoef1()
getWaveletCoef2()
getRidgelineIntensity()
getRidgelineLength()
getCoef2Area()
getPeakShape()