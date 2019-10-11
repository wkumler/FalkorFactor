
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

# Peak object definition ----
peak_object <- setClass("peak_object", slots = list(center="numeric",
                                                    height="numeric",
                                                    width="numeric",
                                                    area="numeric",
                                                    ints="numeric",
                                                    possible_centers = "numeric",
                                                    scan_start="numeric",
                                                    scan_end="numeric",
                                                    gauss_fit="numeric",
                                                    height_top3="numeric",
                                                    ridge_length="numeric",
                                                    ridge_drift="numeric"))



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
  unique(peak_instance@possible_centers[which.max(wavelet_ints)])
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
  xcms:::descendMinTol(roi$int, maxDescOutlier = 2,
                       startpos = c(left_shoulder_offset, right_shoulder_offset))
}

widthFinder <- function(peak_ints, peak_center, peakwidth){
  peak_widths_to_check <- seq(min(peakwidth), max(peakwidth), 2)
  peak_ints_buffered <- c(numeric(max(peakwidth)/2), peak_ints, numeric(max(peakwidth)/2))
  peak_fits <- sapply(peak_widths_to_check, function(pred_peak_width){
    perf_peak <- perf_peak_list[[as.character(pred_peak_width)]]
    peak_left <- (peak_center+max(peakwidth)/2-pred_peak_width/2)
    peak_right <- (peak_center+max(peakwidth)/2+pred_peak_width/2)
    relevant_ints <- peak_ints_buffered[peak_left:peak_right]
    return(cor(relevant_ints, perf_peak))
  })
  peak_width <- peak_widths_to_check[which.max(peak_fits)]
  return(list(edges=c(floor(peak_center-peak_width/2), 
                      ceiling(peak_center+peak_width/2)),
              cor=max(peak_fits)))
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

#' Extend a Region of Interest to account for missed scans
#' 
#' \code{extendROI} accepts a data frame containing ROI information (mz, rt, int)
#' and extends it \code{ext_width} number of seconds in the retention
#' time direction. This accounts for peaks that may be split by missed scans, but
#' also introduces the potential for overlapping peaks and is less likely to
#' throw out peaks because each ROI is more likely to be long enough to 
#' theoretically support a peak.
#' 
#' @param roi_df A Region of Interest data frame, with rt, mz, and int columns.
#' Typically produced by calling \code{split()} on a complete EIC.
#' 
#' @param ext_width The time, in seconds, that the ROI will be extended in both
#' directions. If close to the EIC boundaries, the ROI will not be extended
#' beyond them.
#' 
#' @param eic An Extracted Ion Chromatogram, typically produced by the first
#' portion of the Thanos workflow. A dataframe, with mz, rt, and int information
#' each in its own column.
#' 
#' @return The extended ROI.
extendROI <- function(roi_df, ext_width, eic = eic){
  min_ext_rt <- max(min(eic$rt), min(roi_df$rt)-ext_width)
  max_ext_rt <- min(max(eic$rt), max(roi_df$rt)+ext_width)
  d <- subset(eic, eic$rt>min_ext_rt&eic$rt<max_ext_rt)
}

peakCheck <- function(peak_id){
  idxs <- as.numeric(strsplit(peak_id, "\\.")[[1]])
  peak_data <- all_peak_list[[idxs[1]]][[idxs[2]]][[idxs[3]]]
  plot(peak_data$EIC_rts, peak_data$EIC_ints,
       type="l", lwd=2, 
       xlim=c(min(peak_data$EIC_rts), max(peak_data$EIC_rts)*1.5))
  lines(peak_data$Peak_rts, peak_data$Peak_ints,
        lwd=2, col="red")
  reportvals <- suppressWarnings(sapply(as.numeric(peak_data[sapply(peak_data, length)<=1])[-1], round, digits=2))
  reportnames <- gsub("Peak_", "", names(peak_data[sapply(peak_data, length)<=1])[-1])
  legend("topright", legend = paste0(reportnames, ": ", reportvals), cex = 0.8)
}




# Generate DF, pre-processing ----
mzs <- lapply(x, mz)
mz <- unlist(mzs, use.names = FALSE)
int <- unlist(lapply(x, intensity), use.names = FALSE)
rts <- unname(unlist(lapply(x, rtime)))
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

# Define variables created within loops
min_cwt_length <- 2^(ceiling(log2(max(peakwidth_scans)*6)))
scales <- (floor(peakwidth_scans[1]/2)):ceiling((peakwidth_scans[2]/2))
possible_peakwidths <- min(peakwidth_scans):max(peakwidth_scans)
perf_peak_list <- lapply(possible_peakwidths, function(x){
  exp((-seq(-2, 2, length.out = x+1)^2))
})
names(perf_peak_list) <- as.character(possible_peakwidths)


# Loop over EICs
all_peak_list <- list()
all_peak_ids <- list()
pb <- txtProgressBar(min = 0, max = length(eic_list), style = 3)
for(i in 1:length(eic_list)){
  setTxtProgressBar(pb, i)
  eic <- eic_list[[i]]
  
  # Split EIC by missed scans to turn into ROIs
  roi_encoding <- rle(rts%in%eic$rt)
  peak_lengths <- roi_encoding$lengths[roi_encoding$values]
  roi_all_list <- split(eic, rep(1:length(peak_lengths), times = peak_lengths))
  
  # Remove all ROIs less than a peakwidth long
  roi_list <- roi_all_list[sapply(roi_all_list, nrow)>min(peakwidth_scans)]
  if(!length(roi_list)){ # If there were no reasonable ROIs found
    next
  }
  
  # Loop over ROIs
  rois_per_eic <- list()
  for(j in 1:length(roi_list)){
    roi <- roi_list[[j]]
    roi_start_scan <- which(rts==roi[1, "rt"])-1
    
    # Extend the ROI if necessary to ensure all scales are run
    if(nrow(roi)<min_cwt_length){ #Buffer with zeros at end
      roi_intensity <- c(roi$int, integer(min_cwt_length-nrow(roi)))
    } else {
      roi_intensity <- roi$int
    }
    wcoef_matrix <- xcms:::MSW.cwt(roi_intensity, scales, wavelet = "mexh")
    possible_peaks <- wcoef_matrix[1:nrow(roi),] %>%
      xcms:::MSW.getLocalMaximumCWT() %>%
      xcms:::MSW.getRidge()

    # Loop over peaks
    peaks_per_roi <- list()
    for(k in 1:length(possible_peaks)){
      peak_k <- peak_object(possible_centers = possible_peaks[[k]])
      
      #Center
      peak_k@center <- findPeakCenter(peak_k, roi$int)
      
      #Edges and fit
      edges_output <- widthFinder(roi$int, peak_k@center, peakwidth_scans)
      peak_k@gauss_fit <- edges_output$cor
      peak_k@scan_start <- max(1, edges_output$edges[1])
      peak_k@scan_end <- min(edges_output$edges[2], nrow(roi))
      peak_k@width <- peak_k@scan_end-peak_k@scan_start
      if(peak_k@width < peakwidth_scans[1]){
        next
      }
      
      #Height and intensity
      peak_k@ints <- roi$int[seq(peak_k@scan_start, peak_k@scan_end)]
      peak_k@height <- max(peak_k@ints)
      peak_k@height_top3 <-
        mean(sort(peak_k@ints, partial=peak_k@width-2)[peak_k@width:(peak_k@width-2)])
      peak_k@area <- findPeakArea(peak_k)
      
      #Ridge info
      peak_k@ridge_length <- length(peak_k@possible_centers)
      peak_k@ridge_drift <- length(unique(peak_k@possible_centers))/peak_k@ridge_length
      
      #Signal-to-noise
      #peak_k@SNR <- getSNR(peak_k, eic$int)
      
      
      peaks_per_roi[[k]] <- 
        list("Peak_id"=paste(i, j, k, sep = "."),
              "Peak_center"=rts[peak_k@center+roi_start_scan], 
              "Peak_height"=peak_k@height, 
              "Peak_width"=peak_k@width*mean(diff(rts)), 
              "Peak_area"=peak_k@area, 
              "Peak_start_time"=rts[peak_k@scan_start+roi_start_scan],
              "Peak_end_time"=rts[peak_k@scan_end+roi_start_scan], 
              "Peak_area_top"=peak_k@height_top3, 
              "Peak_ridge_length"=peak_k@ridge_length, 
              "Peak_ridge_drift"=peak_k@ridge_drift,
              "Peak_gauss_fit"=peak_k@gauss_fit,
              "EIC_ints"=eic$int,
              "EIC_rts"=eic$rt,
              "Peak_ints"=peak_k@ints,
              "Peak_rts"=rts[(peak_k@scan_start:peak_k@scan_end)+roi_start_scan])
      all_peak_ids[[length(all_peak_ids)+1]] <- paste(i, j, k, sep = ".")
    }
    rois_per_eic[[j]] <- peaks_per_roi
  }
  all_peak_list[[i]] <- rois_per_eic
}
close(pb)


# Assemble single values into data frame ----
peak_df <- as.data.frame(do.call(rbind, lapply(all_peak_ids, function(x){
  idxs <- as.numeric(strsplit(x, "\\.")[[1]])
  unlist(all_peak_list[[idxs[1]]][[idxs[2]]][[idxs[3]]][sapply(
    all_peak_list[[idxs[1]]][[idxs[2]]][[idxs[3]]], length)<=1])
})))
all_peak_ids <- unlist(all_peak_ids)
# Convert columns to numeric values
for(i in 2:ncol(peak_df)){
  peak_df[,i] <- as.numeric(as.character(peak_df[,i]))
}
peak_df <- arrange(peak_df, desc(Peak_gauss_fit))
#peak_df <- arrange(peak_df, Peak_ridge_drift)

for(i in peak_df$Peak_id){
  peakCheck(i)
  readline(prompt = "Press Enter")
}

