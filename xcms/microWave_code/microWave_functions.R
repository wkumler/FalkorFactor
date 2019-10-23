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
widthFinder <- function(peak_ints, peak_center, peakwidth, perf_peak_list){
  peak_widths_to_check <- seq(min(peakwidth), max(peakwidth), 2)
  peak_ints_buffered <- c(numeric(max(peakwidth)/2), peak_ints, numeric(max(peakwidth)/2))
  peak_fits <- sapply(peak_widths_to_check, function(pred_peak_width, perf_peak_list){
    perf_peak <- perf_peak_list[[as.character(pred_peak_width)]]
    peak_left <- (peak_center+max(peakwidth)/2-pred_peak_width/2)
    peak_right <- (peak_center+max(peakwidth)/2+pred_peak_width/2)
    relevant_ints <- peak_ints_buffered[peak_left:peak_right]
    return(cor(relevant_ints, perf_peak))
  }, perf_peak_list)
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
#' \code{constructEICs} accepts a data frame containing MS information (mz, rt,
#' int) and constructs extracted ion chromatograms (EICs) from the data. To do
#' this, it identifies the highest-intensity data point and collects all data
#' points within the user-specified ppm. All data points that fall within this
#' window are binned into a single EIC and removed from later consideration.
#'
#' The highest-intensity data point is used because accuracy is proportional to
#' mass, and thus represents the most probable EIC center value. This dynamic
#' binning method resolves the problem of setting a single bin width because the
#' created bin widths are a function of compound mass, with larger bins
#' automatically being used with large mass values as proportional to the
#' accuracy of the instrument.
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
constructEICs <- function(given_data_frame, ppm = 2.5, report = TRUE,
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



#' Identify peaks within extracted ion chromatograms
#'
#' \code{microWavePeaks} accepts a list of EICs (columns mz, rt, and int), 
#' typically produced by \code{constructEICs}. It then checks each EIC for
#' potential peaks via xcms's CentWave method and scores them on various
#' metrics, returning a data frame of peaks and their scores, as well
#' as their intensity-weighted m/z values, the centers, widths, and absolute
#' areas.
#' 
#' This function is the core of the peak-identification workflow, and may
#' take a long time to run. Approximately 2,000 peaks were identified in 30
#' seconds from 200 EICs composed of ~100,000 individual data points.
#' 
#' @param eic_list A list of EICs, typically produced by \code{constructEICs}
#' 
#' @param peakwidth A length-two integer vector with the minimum and maximum
#' acceptable peak widths. Defaults to \code{c(20, 80)} for HILIC data
#' 
#' @return A data frame of peak information, with columns
#' \itemize{
#'   \item \strong{Peak_id:} A unique character string to identify each peak,
#'   format "EIC.ROI.Peaknum"
#'   \item \strong{Peak_mz:} The intensity-weighted m/z value of the peak, in Da.
#'   \item \strong{Peak_cener:} The center of the peak, as determined by the 
#'   maximum intensity +/- 5 scans
#'   \item \strong{Peak_height:} The maximum height of the peak. Not necessarily
#'   the intensity at the peak center, as some peak centers are local minima
#'   \item \strong{Peak_width:} The width of the peak, in seconds
#'   \item \strong{Peak_area:} The absolute area of the peak, as calculated
#'   by trapezoidal Riemann sum.
#'   \item \strong{Peak_start_time:} The start time of the peak, in seconds
#'   \item \strong{Peak_end_time:} The end time of the peak, in seconds
#'   \item \strong{Peak_area_top:} The average height of the highest three
#'   values in the peak. If there's strong disagreement between this and the
#'   peak_height value, the peak may be a single spike or otherwise poorly
#'   shaped
#'   \item \strong{Peak_ridge_length:} The length of the "ridge" detected for
#'   the given peak, as determined by xcms' MSW.getRidge(). Longer ridges mean
#'   that many wavelets have a local maxima on a given ridge. However, peaks
#'   that are closely co-eluting may be smoothed into a single peak at high
#'   wavelet scales, and thus only one will have a long ridge while the other
#'   is truncated.
#'   \item \strong{Peak_ridge_drift:} The length of the ridge divided by the
#'   number of scans the ridge is found across. Drift values close to 1 indicate
#'   strong disagreement about the central location of the peak, while drift
#'   values close to zero indicate good agreement as to where the peak is 
#'   located.
#'   \item \strong{Peak_gauss_fit:} The correlation coefficient between the
#'   idealized Gaussian curve used to model the peak and the peak data itself.
#'   This has shown to be a strong predictor of peak quality, with values
#'   close to 1 corresponding to nicely Gaussian peaks and values near zero
#'   indicating essentially random distribution of intensity values
#'   \item \strong{Peak_SNR:} The signal-to-noise ratio of the peak. Calculated
#'   by subtracting the average shoulder value (the outer 3 scans on each side
#'   of the identified peak) from the maximum peak intensity, then dividing
#'   by the standard deviation of the residuals from the Gaussian fit. This
#'   metric essentially produces a z-value for the likelihood that the maximum
#'   peak intensity was drawn from a random sample of noise values.
#' }
microWavePeaks <- function(eic_list, rts, peakwidth = c(20, 80)){
  # Define variables created within loops
  original_data <- do.call(rbind, eic_list)
  time_between_scans <- mean(diff(unique(original_data$rt)))
  peakwidth_scans <- c(floor(peakwidth[1]/time_between_scans), 
                       ceiling(peakwidth[2]/time_between_scans))
  min_cwt_length <- 2^(ceiling(log2(max(peakwidth_scans)*6)))
  scales <- (floor(peakwidth_scans[1]/2)):ceiling((peakwidth_scans[2]/2))
  possible_peakwidths <- min(peakwidth_scans):max(peakwidth_scans)
  perf_peak_list <- lapply(possible_peakwidths, function(x){
    exp((-seq(-2.5, 2.5, length.out = x+1)^2))
  })
  names(perf_peak_list) <- as.character(possible_peakwidths)
  #rts <- unique(original_data$rt)
  
  
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
      wcoef_matrix <- wcoef_matrix[1:nrow(roi),]
      local_maxima <- xcms:::MSW.getLocalMaximumCWT(wcoef_matrix)
      possible_peaks <- xcms:::MSW.getRidge(local_maxima)
      
      # Loop over peaks
      peaks_per_roi <- list()
      for(k in 1:length(possible_peaks)){
        peak_k <- peak_object(possible_centers = possible_peaks[[k]])
        
        #Center
        peak_k@center <- findPeakCenter(peak_k, roi$int)
        
        #Edges, fit, and SNR
        fitting_output <- widthFinder(roi$int, peak_k@center, peakwidth_scans, perf_peak_list)
        peak_k@gauss_fit <- fitting_output$cor
        peak_k@scan_start <- max(1, fitting_output$edges[1])
        peak_k@scan_end <- min(fitting_output$edges[2], nrow(roi))
        peak_k@width <- peak_k@scan_end-peak_k@scan_start
        peak_k@SNR <- fitting_output$SNR
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
        
        #Mz info
        peak_mzs <- roi$mz[seq(peak_k@scan_start, peak_k@scan_end)]
        peak_k@mz <- weighted.mean(peak_mzs, peak_k@ints)
        
        
        peaks_per_roi[[k]] <- 
          list("Peak_id"=paste(i, j, k, sep = "."),
               "Peak_mz"=peak_k@mz,
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
               "Peak_SNR"=peak_k@SNR,
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
  
  peak_df <- as.data.frame(do.call(rbind, lapply(all_peak_ids, function(x){
    idxs <- as.numeric(strsplit(x, "\\.")[[1]])
    unlist(all_peak_list[[idxs[1]]][[idxs[2]]][[idxs[3]]][sapply(
      all_peak_list[[idxs[1]]][[idxs[2]]][[idxs[3]]], length)<=1])
  })))
  for(i in 2:ncol(peak_df)){
    peak_df[,i] <- as.numeric(as.character(peak_df[,i]))
  }
  print(paste("Found", nrow(peak_df), "peaks!"))
  return(peak_df)
}



#' Plot the EIC for a given peak
#' 
#' \code{peakCheck} accepts the output from \code{constructEICs} and
#' \code{microWavePeaks} as well as a specific peak ID, then plots the EIC in
#' which the peak was found with the peak itself highlighted in red. This is
#' often most helpful after sorting the peak data frame by some quality metric
#' then iterating over the top choices to identify where peak quality drops
#' below an arbitrary threshold.
#' 
#' @param eic_list A list of EICs, typically produced by \code{constructEICs}
#' 
#' @param peak_df A data frame of peaks found within the \code{eic_list}. Must
#' include Peak_id, Peak_start_time, Peak_end_time, and Peak_width columns.
#' Typically produced by \code{microWavePeaks}.
#' 
#' @param peak_id A given peak ID, as found in the peak_df data frame. Format
#' should be "EIC.ROI.peaknum".
#' 
#' @param zoom An optional way to zoom in on the peak itself, setting the plot
#' axis limits to frame the peak itself rather than the EIC as a whole. Defaults
#' to FALSE.
#' 
#' @param pts An optional way to include dots at each data point, which sometimes
#' helps to identify missed scans or especially clean peaks that may be missed
#' by other peakfinding software. Defaults to FALSE.
peakCheck <- function(eic_list, peak_df, peak_id, zoom=F, pts=F){
  peak_info <- subset(peak_df, Peak_id==peak_id)
  eic_data <- eic_list[[as.numeric(strsplit(peak_id, "\\.")[[1]])[1]]]
  peak_data <- subset(eic_data, rt>peak_info$Peak_start_time&rt<peak_info$Peak_end_time)
  
  if(zoom){
    ylimits <- c(min(peak_data$int), max(peak_data$int)*1.2)
    xlimits <- c(min(peak_data$rt)-5, min(peak_data$rt)+peak_info$Peak_width*1.5)
  } else {
    xlimits <- c(min(eic_data$rt), max(eic_data$rt)*1.5)
    ylimits <- c(min(eic_data$int), max(eic_data$int)*1.2)
  }
  
  plot(eic_data$rt, eic_data$int, type="l", lwd=2, #pch = 19, cex=0.3,
       xlim=xlimits, ylim=ylimits)
  
  if(pts){
    points(eic_data$rt, eic_data$int, pch=19)
    lines(peak_data$rt, peak_data$int, lwd=2, col="red")
    points(peak_data$rt, peak_data$int, lwd=2, col="red", pch=21, bg="white")
  } else {
    lines(peak_data$rt, peak_data$int, lwd=2, col="red")
  }

  reportvals <- c(peak_info$Peak_mz, suppressWarnings(sapply(as.numeric(peak_info[sapply(peak_info, length)<=1])[-c(1,2)], round, digits=2)))
  reportnames <- gsub("Peak_", "", names(peak_info[sapply(peak_info, length)<=1])[-1])
  legend("topright", legend = paste0(reportnames, ": ", reportvals), cex = 0.8)
  legend("topleft", legend = paste0("Peak id: ", peak_id))
}
