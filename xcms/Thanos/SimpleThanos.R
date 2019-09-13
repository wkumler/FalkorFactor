
# Setup things ----
library(xcms)
library(tidyverse)
library(roxygen2)

load("xcms/raw_data")
x <- raw_data %>%
  filterMsLevel(msLevel. = 1L) %>%
  selectFeatureData(fcol = c(MSnbase:::.MSnExpReqFvarLabels, "centroided")) %>%
  lapply(1:length(fileNames(.)), FUN=filterFile, object = .) %>%
  `[[`(1) %>%
  spectra()

# Peak objects ----
peak_object <- setClass("peak_object", slots = list(center="numeric",
                                                    height="numeric",
                                                    width="numeric",
                                                    area="numeric",
                                                    ints="numeric",
                                                    wavelet_coefs="matrix",
                                                    local_maxima="matrix",
                                                    possible_centers = "numeric",
                                                    num_used_scales="numeric",
                                                    best_scale="numeric",
                                                    scan_start="numeric",
                                                    scan_end="numeric",
                                                    height_top3="numeric",
                                                    coef2area="numeric",
                                                    ridge_length="numeric",
                                                    ridge_prop="numeric",
                                                    ridge_drift="numeric",
                                                    linearity="numeric",
                                                    coef_fit="numeric",
                                                    sigma_star="numeric",
                                                    norm_sigma_star="numeric"))



# Peak functions ----

#' Finds local maxima by summing nearby intensities
#' 
#' \code{findPeakCenter} accepts a peak object and the intensities of the ROI
#' containing the peak object. The peak object must have a filled "possible
#' centers" slot containing a vector of scans, usually produced by calling
#' the xcms:::MSW.getRidge() function on the peak object's local max matrix.
#' The function then iterates over the possible centers, calculates the
#' area of 5 scans around the center, and returns the scan with the highest
#' area.
#' 
#' @param peak_instance The peak_instance parameter is an S4 peak object, 
#' usually produced by the normal Thanos workflow.
#' 
#' @param roi_ints The roi_ints parameter is a vector of intensities 
#' corresponding to the ROI. These intensities will be used to calculate the
#' area around each possible peak.
#' 
#' @return The scan number corresponding to the highest area - i.e., the center
#' of the peak.
findPeakCenter <- function(peak_instance, roi_ints){
  if(!length(peak_instance@possible_centers)){ # If slot is empty
    stop("The 'possible_centers' slot of this object is empty")
  }
  wavelet_ints <- sapply(unique(peak_instance@possible_centers), function(x){
    sum(roi_ints[max(1, x-5):(x+5)], na.rm = T)
  })
  unique(peak_instance@possible_centers)[which.max(wavelet_ints)]
}

#' Finds the best wavelet scale for a given peak object
#' 
#' \code{findBestScale} accepts a peak object which must have filled 
#' "local_maxima" and "center" slots. The function then finds the largest
#' wavelet scale which has a local maximum in the same scan as the peak center.
#' This function should under go serious revision, as this algorithm's utility
#' is questionable at best.
#' 
#' @param peak_instance The peak_instance parameter is an S4 peak object, 
#' usually produced by the normal Thanos workflow. Slots "local_maxima" and
#' "center" must be filled.
#' 
#' @return The wavelet scale corresponding to the wavelet sharing a local
#' maximum with the peak center
findBestScale <- function(peak_instance){
  if(!length(peak_instance@local_maxima)){ # If slot is empty
    stop("The 'local_maxima' slot of this object is empty")
  }
  if(!length(peak_instance@center)){ # If slot is empty
    stop("The 'center' slot of this object is empty")
  }
  scale_maxes <- as.logical(peak_instance@local_maxima[peak_instance@center,])
  max(as.numeric(colnames(local_maxima))[scale_maxes])
}

#' Finds the edges of a given peak object
#' 
#' \code{findPeakEdges} accepts a peak object which must have filled 
#' "best_scale" and "center" slots. The function then finds the edges of the
#' peak by using xcms:::descendMinTol, with center +/- best_scale as the
#' offset parameters.
#' 
#' @param peak_instance The peak_instance parameter is an S4 peak object, 
#' usually produced by the normal Thanos workflow. Slots "best_scale" and
#' "center" must be filled.
#' 
#' @return A vector with two values, corresponding to the left and right edges
#' of the peak, respectively.
findPeakEdges <- function(peak_instance){
  left_shoulder_offset <- peak_instance@center-peak_instance@best_scale
  right_shoulder_offset <- peak_instance@center+peak_instance@best_scale
  xcms:::descendMinTol(roi$int, maxDescOutlier = 0,
                       startpos = c(left_shoulder_offset, right_shoulder_offset))
}

#' Find the total area underneath a curve by Riemann trapezoidal sum
#' 
#' \code{findPeakArea} accepts a peak object with pre-established start, end,
#' width, and intensities. It then calculates the area underneath the peak
#' between the scan's start and end using an exact Riemann trapezoidal sum
#' method.
#' 
#' @param peak_instance The peak_instance parameter is an S4 peak object, 
#' usually produced by the normal Thanos workflow. Slots "scan_start",
#' "peak_width", "peak_ints", and "scan_end" must be filled.
#' 
#' @return A single number, the area under the curve.
findPeakArea <- function(peak_instance){
  given_peak_scans <- peak_instance@scan_start:peak_instance@scan_end
  idx = 2:peak_instance@width
  as.numeric(((given_peak_scans[idx] - given_peak_scans[idx-1]) %*% 
                (peak_instance@ints[idx] + peak_instance@ints[idx-1]))/2)
}


# Generate DF, pre-processing ----
mzs <- lapply(x, mz)
mz <- unlist(mzs, use.names = FALSE)
int <- unlist(lapply(x, intensity), use.names = FALSE)
rts <- unlist(lapply(x, rtime))
rt <- rep(rts, sapply(mzs, length))
all_data <- data.frame(mz, int, rt)



# Find all reasonable EICs ----
peakwidth <- c(20, 80)
peakwidth_scans <- c(floor(peakwidth[1]/mean(diff(rts))), 
                     ceiling(peakwidth[2]/mean(diff(rts))))
prefilter <- c(1, 1)
ppm <- 2.5

data <- all_data %>% filter(mz>100&mz<120) %>% filter(rt>60&rt<1100)

pb <- txtProgressBar(min = 0, max = 1, style = 3)
eic_list <- list()
data_start_length <- nrow(data)
while(nrow(data)>0){
  setTxtProgressBar(pb, 1-(sqrt(nrow(data))/sqrt(data_start_length)))
  point_of_interest <- data[which.max(data$int),]
  epsilon_Da <- point_of_interest$mz*ppm/1000000
  upper_eic_mz <- point_of_interest$mz+epsilon_Da
  lower_eic_mz <- point_of_interest$mz-epsilon_Da
  eic <- filter(data, mz>lower_eic_mz & mz<upper_eic_mz)
  data <- data[data$mz<lower_eic_mz | data$mz>upper_eic_mz,]
  
  # If the EIC can't contain a peak bc too short, remove it
  if(nrow(eic)<peakwidth_scans[1]){
    next
  }
  
  # If there aren't enough points above the prefilter threshold
  runs_above_prefilter <- rle(eic$int>prefilter[2])
  max_run_length <- max(runs_above_prefilter$lengths[runs_above_prefilter$values])
  if(max_run_length<prefilter[1]){
    next
  }
  
  # Put it all together and save
  eic_list[[length(eic_list)+1]] <- eic
}
close(pb)
print(paste("Extracted", length(eic_list), "ion chromatograms! Now finding peaks..."))



# Find peaks ----
# Loop over EICs
all_peak_list <- list()
pb <- txtProgressBar(min = 0, max = length(eic_list), style = 3)
for(i in 1:length(eic_list)){
  setTxtProgressBar(pb, i)
  eic <- eic_list[[i]]
  # Split EIC by missed scans to turn into ROIs
  roi_encoding <- rle(rts%in%eic$rt)
  peak_lengths <- roi_encoding$lengths[roi_encoding$values]
  roi_all_list <- split(eic, rep(1:length(peak_lengths), times = peak_lengths))
  roi_list <- roi_all_list[sapply(roi_all_list, nrow)>peakwidth_scans[1]]
  if(!length(roi_list)){ # If there were no reasonable ROIs found
    next
  }
  
  # Loop over ROIs
  eic_peak_list <- list()
  for(j in 1:length(roi_list)){
    roi <- roi_list[[j]]
    roi_start_scan <- which(rts==roi[1, "rt"])-1
    
    # Make some waves (same for all peaks in ROI)
    scales <- 1:(peakwidth_scans[2]/2)
    wcoef_matrix <- xcms:::MSW.cwt(roi$int, scales, wavelet = "mexh")
    if(length(wcoef_matrix)==1){ # If CWT returns NA because the scales suck
      next
    }
    local_maxima <- xcms:::MSW.getLocalMaximumCWT(wcoef_matrix)
    possible_peaks <- xcms:::MSW.getRidge(local_maxima)
    num_scales <- length(attr(possible_peaks, "scales"))
    # Remove all peaks that have maxima in less than half the scales ADJUST LATER
    possible_peaks <- possible_peaks[sapply(possible_peaks, length)>(num_scales/2)]
    
    # Loop over peaks
    roi_peak_list <- list()
    for(k in 1:length(possible_peaks)){
      peak_k <- peak_object(wavelet_coefs = wcoef_matrix, 
                            local_maxima = local_maxima,
                            possible_centers = possible_peaks[[k]])
      
      peak_k@center <- findPeakCenter(peak_k, roi$int)
      
      peak_k@best_scale <- findBestScale(peak_k)
      peak_k@num_used_scales <- num_scales
      
      peak_edges <- findPeakEdges(peak_k)
      peak_k@scan_start <- max(1, peak_edges[1])
      peak_k@scan_end <- min(peak_edges[2], nrow(roi))
      peak_k@width <- peak_k@scan_end-peak_k@scan_start
      
      if(peak_k@width < peakwidth_scans[1] || peak_k@width > peakwidth[2]){
        next
      }
      
      peak_ints <- roi$int[seq(peak_k@scan_start:peak_k@scan_end)]
      peak_k@height <- max(peak_ints)
      
      # peak_k@height_top3 <- 
      #   mean(sort(peak_ints, partial=peak_k@width-2)[peak_k@width:(peak_k@width-2)])
      
      peak_k@area <- findPeakArea(peak_k)
      
      peak_coefs <- peak_k@wavelet_coefs[peak_k@scan_start:peak_k@scan_end, ]
      peak_k@coef2area <- max(peak_coefs)/peak_k@area
      
      peak_k@ridge_length <- length(peak_k@possible_centers)
      peak_k@ridge_prop <- round(peak_k@ridge_length/peak_k@num_used_scales, digits = 2)
      peak_k@ridge_drift <- length(unique(peak_k@possible_centers))/peak_k@ridge_length
      
      peak_k@linearity <- cor(sort(peak_ints), 1:length(peak_ints))^2
      
      peak_k@coef_fit <- cor(peak_ints, peak_coefs[, peak_k@best_scale])
      
      peak_mzs <- roi$mz[peak_k@scan_start:peak_k@scan_end]
      peak_k@sigma_star <- max(peak_mzs)-min(peak_mzs)
      peak_k@norm_sigma_star <- peak_k@sigma_star/(max(roi$mz)-min(roi$mz))
      
      roi_peak_list[[k]] <- 
        cbind("Peak_id"=paste(i, j, k, sep = "."),
              "Peak_center"=rts[peak_k@center+roi_start_scan], 
              "Peak_height"=peak_k@height, 
              "Peak_width"=peak_k@width*mean(diff(rts)), 
              "Peak_area"=peak_k@area, 
              "Peak_best_scale"=peak_k@best_scale, 
              "Peak_start_time"=rts[peak_k@scan_start+roi_start_scan],
              "Peak_end_time"=rts[peak_k@scan_end+roi_start_scan], 
              #"Peak_area_top"=peak_k@height_top3, 
              "Peak_coef2area"=peak_k@coef2area, 
              "Peak_ridge_length"=peak_k@ridge_length, 
              "Peak_ridge_prop"=peak_k@ridge_prop, 
              "Peak_ridge_drift"=peak_k@ridge_drift,
              "Peak_linearity"=peak_k@linearity,
              "Peak_coef_fit"=peak_k@coef_fit,
              "Peak_sigma_star"=peak_k@sigma_star,
              "Peak_norm_sigma_star"=peak_k@norm_sigma_star)
    }
    roi_peak_df <- as.data.frame(do.call(rbind, roi_peak_list))
    eic_peak_list[[j]] <- roi_peak_df
  }
  eic_peak_df <- do.call(rbind, eic_peak_list)
  all_peak_list[[i]] <- eic_peak_df
}
close(pb)
all_peak_df <- do.call(rbind, all_peak_list)

