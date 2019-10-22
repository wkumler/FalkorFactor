# MicroWave functions

library(roxygen2)

#' An S4 class for holding peak information.
#'
#' @slot center A length-one numeric vector with the retention time  (in
#'   seconds) of a given peak, as determined by the scan with the highest
#'   intensity +/- 5 scans
#'
#' @slot height A length-one numeric vector with the maximum height of the peak
#'
#' @slot mz A length-one numeric vector with the mass-to-charge (m/z) value of
#'   the peak, as weighted by the intensity of each data point
#'
#' @slot width A length-one numeric vector with the width of the peak, in
#'   seconds
#'
#' @slot area A length-one numeric vector with the absolute area of the peak, as
#'   determined by the exact Riemann sum of intensities
#'
#' @slot ints A numeric vector with the intensity values of a peak
#'
#' @slot possible_centers A numeric vector with the possible centers of a peak,
#'   one of which will later become the actual center after processing with
#'   findPeakCenter
#'
#' @slot scan_start A length-one numeric vector with the start scan number of
#'   the peak
#'
#' @slot scan_end A length-one numeric vector with the final scan number of the
#'   peak
#'
#' @slot gauss_fit A length-one numeric vector between 0 and 1 with the
#'   correlation coefficient between the peak intensities and the idealized
#'   Gaussian curve with a width equal to the width of the peak
#'
#' @slot height_top3 A length-one numeric vector with the average intensity of
#'   the three highest intensities of the peak
#'
#' @slot ridge_length A length-one numeric vector with the length of the ridge
#'   corresponding to the given peak, as determined by getRidge()
#'
#' @slot ridge_drift A length-one numeric vector between 0 and 1 with the
#'   "drift" of the ridge. This is calculated by dividing the number of unique
#'   ridge maxima (in scans) by the length of the ridge. Higher drifts imply
#'   disagreement between the wavelets about the location of the ridge maximum,
#'   and lower drifts indicate that the ridge maxima all occur at the same scan.
#'
#' @slot SNR A length-one numeric vector with the signal-to-noise ratio of the
#'   peak. This is calculated by dividing the maximum peak intensity minus the
#'   baseline intensity by the standard deviation of the noise within the peak.
#'   That peak noise is calculated from the residuals of the Gaussian
#'   correlation, measuring the deviation from the ideal curve. This assumes
#'   that the peak is Gaussian and that the noise within the peak is the same as
#'   the noise outside the peak after the peak signal has been subtracted.
peak_object <- setClass("peak_object", slots = list(center="numeric",
                                                    height="numeric",
                                                    mz="numeric",
                                                    width="numeric",
                                                    area="numeric",
                                                    ints="numeric",
                                                    possible_centers = "numeric",
                                                    scan_start="numeric",
                                                    scan_end="numeric",
                                                    gauss_fit="numeric",
                                                    height_top3="numeric",
                                                    ridge_length="numeric",
                                                    ridge_drift="numeric",
                                                    SNR="numeric"))



#' Finds local maxima by summing nearby intensities
#'
#' \code{findPeakCenter} accepts a peak object and the intensities of the ROI
#' containing the peak object. The peak object must have a filled "possible
#' centers" slot containing a vector of scans, usually produced by calling the
#' xcms:::MSW.getRidge() function on the peak object's local max matrix. The
#' function then iterates over the possible centers, calculates the area of 5
#' scans around the center, and returns the scan with the highest area.
#'
#' @param peak_instance The peak_instance parameter is an S4 peak object,
#'   usually produced by the normal Thanos workflow.
#'
#' @param roi_ints The roi_ints parameter is a vector of intensities
#'   corresponding to the ROI. These intensities will be used to calculate the
#'   area around each possible peak.
#'
#' @return The scan number corresponding to the highest area - i.e., the center
#'   of the peak.
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
#'   usually produced by the normal Thanos workflow. Slots "local_maxima" and
#'   "center" must be filled.
#'
#' @return The wavelet scale corresponding to the wavelet sharing a local
#'   maximum with the peak center
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
#' peak by using xcms:::descendMinTol, with center +/- best_scale as the offset
#' parameters. This code is not currently being used by the microWave workflow,
#' as it's been replaced by widthFinder instead.
#'
#' @param peak_instance The peak_instance parameter is an S4 peak object,
#'   usually produced by the normal Thanos workflow. Slots "best_scale" and
#'   "center" must be filled.
#'
#' @return A vector with two values, corresponding to the left and right edges
#'   of the peak, respectively.
findPeakEdges <- function(peak_instance){
  left_shoulder_offset <- peak_instance@center-peak_instance@best_scale
  right_shoulder_offset <- peak_instance@center+peak_instance@best_scale
  xcms:::descendMinTol(roi$int, maxDescOutlier = 2,
                       startpos = c(left_shoulder_offset, right_shoulder_offset))
}


#' Finds the edges of a peak, fits a Gaussian curve to it, and calculates SNR
#'
#' \code{widthFinder} accepts peak information (intensity, peak center, and
#' maximum/minimum peakwidths) and uses a Gaussian fitting model to calculate
#' the edges of the peak. In essence, the model fits a Gaussian curve
#' corresponding to each possible peakwidth and tests which curve fits the peak
#' best. 3 SDs to the left and right are used as peak edges, as long as they're
#' within a contiguous set of scans. The Gaussian fit of the best curve is
#' reported, which appears to be a robust parameter for peak identification.
#' Finally, the residuals of the Gaussian curve vs the actual peak are
#' calculated and used to estimate the noise (curve - normalized peak
#' intensity)/SD(residuals) and the three outermost scans on the left and the
#' right are used to estimate the background intensity level.
#'
#' @param peak_ints A vector of intensities corresponding to a peak. These must
#'   be contiguous (no missed scans) and typically correspond to all the
#'   intensities within an ROI.
#'
#' @param peak_center A single value, obtained from the wavelet matching
#'   algorithm in XCMS.
#'
#' @param peakwidth_scans The possible peak widths in scans. This value is
#'   typically set early on in the workflow in terms of seconds and converted to
#'   scans using the average time between scans.
#'
#' @return A list containing the peak edges, the correlation coefficient of the
#'   Gaussian fit, and the signal-to-noise parameter.
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
  best_peak_width <- peak_widths_to_check[which.max(peak_fits)]
  
  peak_left <- (peak_center+max(peakwidth)/2-best_peak_width/2)
  peak_right <- (peak_center+max(peakwidth)/2+best_peak_width/2)
  relevant_ints <- peak_ints_buffered[peak_left:peak_right]
  
  best_perf_peak <- perf_peak_list[[best_peak_width-min(peakwidth)+1]]
  
  residuals <- best_perf_peak-relevant_ints/max(relevant_ints)
  peak_noise <- sd(residuals)*max(peak_ints)
  peak_background <- mean(c(head(peak_ints, 3), tail(peak_ints, 3)))
  SNR <- max(peak_ints)/(peak_background+peak_noise)
  
  return(list(edges=c(floor(peak_center-best_peak_width/2), 
                      ceiling(peak_center+best_peak_width/2)),
              cor=max(peak_fits),
              SNR=SNR))
}

#' Find the total area underneath a curve by Riemann trapezoidal sum
#'
#' \code{findPeakArea} accepts a peak object with pre-established start, end,
#' width, and intensities. It then calculates the area underneath the peak
#' between the scan's start and end using an exact Riemann trapezoidal sum
#' method.
#'
#' @param peak_instance The peak_instance parameter is an S4 peak object,
#'   usually produced by the normal Thanos workflow. Slots "scan_start",
#'   "peak_width", "peak_ints", and "scan_end" must be filled.
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
#' \code{extendROI} accepts a data frame containing ROI information (mz, rt,
#' int) and extends it \code{ext_width} number of seconds in the retention time
#' direction. This accounts for peaks that may be split by missed scans, but
#' also introduces the potential for overlapping peaks and is less likely to
#' throw out peaks because each ROI is more likely to be long enough to
#' theoretically support a peak.
#'
#' @param roi_df A Region of Interest data frame, with rt, mz, and int columns.
#'   Typically produced by calling \code{split()} on a complete EIC.
#'
#' @param ext_width The time, in seconds, that the ROI will be extended in both
#'   directions. If close to the EIC boundaries, the ROI will not be extended
#'   beyond them.
#'
#' @param eic An Extracted Ion Chromatogram, typically produced by the first
#'   portion of the Thanos workflow. A dataframe, with mz, rt, and int
#'   information each in its own column.
#'
#' @return The extended ROI.
extendROI <- function(roi_df, ext_width, eic = eic){
  min_ext_rt <- max(min(eic$rt), min(roi_df$rt)-ext_width)
  max_ext_rt <- min(max(eic$rt), max(roi_df$rt)+ext_width)
  d <- subset(eic, eic$rt>min_ext_rt&eic$rt<max_ext_rt)
}


#' Construct extracted ion chromatograms from raw data values
#'
#' \code{makeEICs} accepts a data frame containing MS information (mz, rt, int)
#' and constructs extracted ion chromatograms (EICs) from the data. To do this,
#' it identifies the highest-intensity data point and collects all data points
#' within the user-specified ppm. All data points that fall within this window
#' are binned into a single EIC and removed from later consideration.
#'
#' The highest-intensity data point is used because accuracy is inversely
#' proportional to mass, and thus represents the most probable EIC center value.
#' This dynamic binning method resolves the problem of setting a single bin
#' width because the created bin widths are a function of compound mass, with
#' larger bins automatically being used with large mass values as proportional
#' to the accuracy of the instrument.
#'
#' @param given_data_frame A raw data frame with rt, mz, and int named columns.
#'
#' @param ppm The parts-per-million accuracy of the instrument on which the data
#'   was collected. Current default is 2.5 for Orbitrap data.
#'
#' @param report A logical value telling the function whether or not to report
#'   the initial number of data points, the progress of the function as
#'   displayed by a text-based loading bar, and the final number of EICs
#'   constructed by the function. The loading bar is set to report the *square
#'   root* of the number of remaining raw data points to provide a more
#'   realistic time estimate, as many points are removed each iteration early on
#'   in the function as EICs which contain data from every scan, while later
#'   iterations take much longer because the EICs are fragmented. Should be set
#'   to FALSE if running in parallel, and the parallel progress bar should be
#'   used instead.
#'
#' @param prefilter A named length-two vector containing prefilter used to
#'   remove EICs that are either too short or fail to have a sufficient number
#'   of continuous data points above the user-specified intensity. The
#'   "contiguous" value is the number of consecutive scans, and the "intensity"
#'   value is the intensity threshold which all contiguous scans must exceed.
#'   Judicious use of the prefilter can significantly speed up the EIC
#'   construction algorithm, but risks removing true peaks if enabled and set
#'   too aggressively. Its functionality can be disabled by setting both values
#'   to zero, as is the default.
#'
#' @param peakwidth A length-two vector containing the minimum and maximum
#'   possible peak width retention times, in seconds. These values are later
#'   converted to minimum and maximum scans by multiplying by the average time
#'   between scans. Defaults to c(20, 80) for HILIC data.
#'
#' @return A list of data frames, each one containing mz, int, and rt values
#'   corresponding to a given EIC.
makeEICs <- function(given_data_frame, ppm = 2.5, report = TRUE,
                     prefilter = c(intensity=0, contiguous=0), 
                     peakwidth = c(20, 80)){
  if(report){print(paste("Constructing EICs from", 
                         nrow(given_data_frame), "data points"))}
  if(report){pb <- txtProgressBar(min = 0, max = 1, style = 3)}
  
  time_between_scans <- mean(diff(unique(given_data_frame$rt)))
  
  peakwidth_scans <- c(floor(peakwidth[1]/time_between_scans), 
                             ceiling(peakwidth[2]/time_between_scans))
  
  eic_list <- list()
  data_start_length <- nrow(given_data_frame)
  while(nrow(given_data_frame)>0){
    if(report){setTxtProgressBar(pb, 1-(sqrt(nrow(given_data_frame))/
                                          sqrt(data_start_length)))}
    point_of_interest <- given_data_frame[which.max(given_data_frame$int),]
    epsilon_Da <- point_of_interest$mz*ppm/1000000
    upper_eic_mz <- point_of_interest$mz+epsilon_Da
    lower_eic_mz <- point_of_interest$mz-epsilon_Da
    eic <- filter(given_data_frame, mz>lower_eic_mz & mz<upper_eic_mz)
    given_data_frame <- given_data_frame[given_data_frame$mz<lower_eic_mz | 
                                           given_data_frame$mz>upper_eic_mz,]
    
    # If the EIC can't contain a peak bc too short, remove it
    if(nrow(eic)<peakwidth_scans[1]){
      next
    }
    
    # If there aren't enough points above the prefilter threshold
    if(any(as.logical(prefilter))){
      runs_above_prefilter <- rle(eic$int>prefilter["intensity"])
      max_run_length <- max(runs_above_prefilter$lengths[runs_above_prefilter$values])
      if(max_run_length<prefilter["contiguous"]){
        next
      }
    }
    
    # Put it all together and save
    eic_list[[length(eic_list)+1]] <- eic
  }
  if(report){close(pb)}
  if(report){print(paste("Extracted", length(eic_list), "ion chromatograms!"))}
  
  return(eic_list)
}
