library(xcms)
library(tidyverse)

load("xcms/raw_data")

object <- raw_data

object_mslevel <- filterMsLevel(
  selectFeatureData(object, fcol = c(MSnbase:::.MSnExpReqFvarLabels, 
                                     "centroided")), msLevel. = 1L)

object_mslevel <- lapply(1:length(fileNames(object_mslevel)), 
                         FUN = filterFile, object = object_mslevel)
object <- object_mslevel[[1]]

x = spectra(object)

mzs <- lapply(x, mz)
mz <- unlist(mzs, use.names = FALSE)

int <- unlist(lapply(x, intensity), use.names = FALSE)

rt <- unlist(lapply(x, rtime), use.names = FALSE)
scantime = rt

vals_per_spect <- lengths(mzs, FALSE)
valsPerSpect = vals_per_spect


peakwidth <- c(20, 80) # Minimum and maximum width of a few random peaks, in seconds
mz_span <- c(0.0005) # Maximum spread of m/z values across a well-defined peak, plus some buffer
ppm <- ceiling((mz_span*1000000)/132)
snthresh = 10
prefilter = c(3, 100000)
mzCenterFun = "weighted.mean"
integrate = 1
mzdiff = -0.001
fitgauss = FALSE
noise = 0
verboseColumns = FALSE
roiList = list()
roiScales = NULL
sleep = 0





# .milliWave2 ----
# The second one
# .milliWave <- function (mz, int, scantime, valsPerSpect, ppm = 25,
#                         peakwidth = c(20, 50), snthresh = 10,
#                         prefilter = c(3, 100), mzCenterFun = "weighted.mean",
#                         integrate = 1, mzdiff = -0.001, fitgauss = FALSE,
#                         noise = 0, sleep = 0, verboseColumns = FALSE,
#                         roiList = list(), firstBaselineCheck = TRUE,
#                         roiScales = NULL)
# {
  valCount <- cumsum(valsPerSpect)
  scanindex <- as.integer(c(0, valCount[-length(valCount)]))
  basenames <- c("mz", "mzmin", "mzmax", 
                 "rt", "rtmin", "rtmax", "into", 
                 "intb", "maxo", "xcms_sn")
  verbosenames <- c("egauss", "mu", "sigma", 
                    "h", "f", "dppm", "scale", "scpos", 
                    "scmin", "scmax", "lmin", "lmax")
  if(length(peakwidth)!=2||any(peakwidth==0)){
    stop("'Peakwidth' must be a non-zero vector of length 2")
  }
  scalerange <- round((peakwidth/mean(diff(scantime)))/2)
  scales <- scalerange[1]:scalerange[2]
  minPeakWidth <- min(scales)
  noiserange <- c(min(scales) * 3, max(scales) * 3)
  minCentroids <- minPtsAboveBaseLine <- max(4, minPeakWidth - 2)
  scRangeTol <- maxDescOutlier <- floor(minPeakWidth/2)
  scanrange <- c(1, length(scantime))
  Nscantime <- length(scantime)
  
  # Find ROIs
  message("Detecting mass traces at ", ppm, " ppm ... ", appendLF = FALSE)
  roiList <- .Call("findmzROI", 
                   mz, int, scanindex, 
                   as.double(c(0, 0)), 
                   as.integer(scanrange), 
                   as.integer(length(scantime)), 
                   as.double(ppm * 1e-06), 
                   as.integer(minCentroids), 
                   as.integer(prefilter), 
                   as.integer(noise), PACKAGE = "xcms")
  message("OK")
  
  
  
  
  
  # Detect peaks
  peaklist <- list()
  message("Detecting chromatographic peaks in ", length(roiList), 
          " regions of interest ...", appendLF = FALSE)
  pb <- txtProgressBar(min = 0, max = length(roiList), style = 3)
  
  #for (f in 1:length(roiList)) { # For each region of interest
    # setTxtProgressBar(pb, f)
  f <- 500
    # Generate EICs
    feat <- roiList[[f]]
    peaks <- peakinfo <- NULL
    mzrange <- c(feat$mzmin, feat$mzmax)
    sccenter <- feat$scmin[1] + round((feat$scmax - feat$scmin)/2)
    scrange <- c(feat$scmin, feat$scmax)
    sr <- c(max(scanrange[1], feat$scmin - max(noiserange)), 
            min(scanrange[2], feat$scmax + max(noiserange)))
    eic <- .Call("getEIC", mz, int, scanindex, as.double(mzrange), 
                 as.integer(sr), as.integer(length(scanindex)), PACKAGE = "xcms")

    d <- eic$intensity
    td <- eic$scan
    
    idxs <- which(eic$scan %in% seq(scrange[1], scrange[2]))
    mzROI.EIC <- list(scan = eic$scan[idxs], intensity = eic$intensity[idxs])
    omz <- .Call("getMZ", mz, int, scanindex, as.double(mzrange), 
                 as.integer(scrange), as.integer(length(scantime)), 
                 PACKAGE = "xcms")
    od <- mzROI.EIC$intensity
    otd <- mzROI.EIC$scan
    if (all(mzROI.EIC$intensity == 0)) {
      warning("centWave: no peaks found in ROI.")
      next
    }
    
    ftd <- max(td[1], scrange[1] - scRangeTol):min(td[length(td)], scrange[2] + scRangeTol)
    fd <- d[match(ftd, td)]
    
    # plot_max <- max(d)
    # plot(td, d, ylim=c(0, plot_max*1.5))
    # arrows(x0 = min(td), x1 = max(td), y0=plot_max*1.03, y1=plot_max*1.03,
    #        code = 3, col = "black", angle = 90)
    # arrows(x0 = min(ftd), x1 = max(ftd), y0=plot_max*1.05, y1=plot_max*1.05,
    #        code = 3, col = "blue", angle = 90)
    # arrows(x0 = min(otd), x1 = max(otd), y0=plot_max*1.07, y1=plot_max*1.07,
    #        code = 3, col = "red", angle = 90)
    # legend("top", legend = c("ROI", "ROI +/- min(peakwidth)/2", "ROI +/- noise"),
    #        col = c("red", "blue", "black"), lwd = 2)
    
    
    
    if ((feat$scmax - feat$scmin + 1) >= 10 * minPeakWidth) {
      noised <- .Call("getEIC", mz, int, scanindex, 
                      as.double(mzrange), as.integer(scanrange), 
                      as.integer(length(scanindex)), 
                      PACKAGE = "xcms")$intensity
    } else {
      noised <- d
    }
    noise <- xcms:::estimateChromNoise(noised, trim = 0.05, minPts = 3 * minPeakWidth)
    nonzero_scans <- noised[noised>0]
    minimum_roi_length <- 3*minPeakWidth
    if(length(nonzero_scans) < minimum_roi_length){ # If roi is dangerously short, avoid trimming
      noise <- mean(noised) # calculate mean INCLUDING zeros?????
    } else {
      noise <- mean(nonzero_scans, trim = 0.05) # Include some, but not all zeros?
    }
    
    
    # Filter #1 below
    roi_rle <- rle(fd>noise)
    max_run_length <- max(roi_rle$lengths[roi_rle$values])
    
    # if (!xcms:::continuousPtsAboveThreshold(fd, threshold = noise, 
    #                                         num = minPtsAboveBaseLine) & 
    #     firstBaselineCheck) { # If the peak doesn't have minPeakWidth-2 scans above average noise
    # if(max_run_length<minPtsAboveBaseLine & firstBaselineCheck){
    #   next # Discard ROI
    # }
    
    lnoise <- xcms:::getLocalNoiseEstimate(d, td, ftd, noiserange, 
                                           length(scantime), threshold = noise, 
                                           num = minPtsAboveBaseLine)
    baseline <- max(1, min(lnoise[1], noise))
    sdnoise <- max(1, lnoise[2])
    sdthr <- sdnoise * snthresh
    # Filter #2
    # if (!(any(fd - baseline >= sdthr))) 
    #   next
    
    wCoefs <- xcms:::MSW.cwt(d, scales = scales, wavelet = "mexh")
    # Filter #3
    # if (!(!is.null(dim(wCoefs)) && any(wCoefs - baseline >= sdthr))) {
    #   next
    # }
    if (length(td) == length(scantime))
      wCoefs[nrow(wCoefs), ] <- wCoefs[nrow(wCoefs) - 1, ] * 0.99
    localMax <- xcms:::MSW.getLocalMaximumCWT(wCoefs)
    rL <- xcms:::MSW.getRidge(localMax)
    
    # Filter #4 (if it doesn't pass this, fails to get passed to the rest of the algorithm)
    wpeaks <- sapply(rL, function(x) {
      w <- min(1:length(x), ncol(wCoefs))
      any(wCoefs[x, w] - baseline >= sdthr)
    })
    
    
    
    if (any(wpeaks)) {
      wpeaksidx <- which(wpeaks)
      for (p in 1:length(wpeaksidx)) { # For each ridgeline detected in an ROI
        opp <- rL[[wpeaksidx[p]]]
        pp <- unique(opp)
        if (length(pp) >= 1) {
          dv <- td[pp] %in% ftd
          if (any(dv)) {
            if (any(d[pp[dv]] - baseline >= sdthr)) { # Filter #5 ?
              if (length(roiScales) > 0) { # Set at the very beginning
                best.scale.nr <- which(scales == roiScales[[f]])
                if (best.scale.nr > length(opp)) 
                  best.scale.nr <- length(opp)
              } else {
                inti <- numeric(length(opp))
                irange <- ceiling(scales[1]/2)
                for (k in 1:length(opp)) {
                  kpos <- opp[k]
                  r1 <- ifelse(kpos - irange > 1, kpos - irange, 1)
                  r2 <- ifelse(kpos + irange < length(d), kpos + irange, length(d))
                  inti[k] <- sum(d[r1:r2])
                }
                best.scale.nr <- which.max(inti)
              }
              best.scale <- scales[best.scale.nr]
              best.scale.pos <- opp[best.scale.nr]
              pprange <- min(pp):max(pp)
              lwpos <- max(1, best.scale.pos - best.scale)
              rwpos <- min(best.scale.pos + best.scale, length(td))
              p1 <- match(td[lwpos], otd)[1]
              p2 <- match(td[rwpos], otd)
              p2 <- p2[length(p2)]
              if (is.na(p1)) 
                p1 <- 1
              if (is.na(p2)) 
                p2 <- feat$scmax - feat$scmin + 1
              mz.value <- omz[p1:p2]
              mz.int <- od[p1:p2]
              maxint <- max(mz.int)
              mzrange <- range(mz.value)
              mzmean <- do.call(mzCenterFun, list(x = mz.value, w = mz.int))
              peaks <- rbind(peaks, c(mzmean, mzrange, 
                                      NA, NA, NA, NA, NA, 
                                      maxint, round((maxint - baseline)/sdnoise), 
                                      NA, NA, NA, NA, 
                                      f, NA, best.scale, td[best.scale.pos], 
                                      td[lwpos], td[rwpos], NA, NA))
              peakinfo <- rbind(peakinfo, c(best.scale, 
                                            best.scale.nr, best.scale.pos, lwpos, 
                                            rwpos))
            }
          }
        }
      }
    }
    
    
    if (!is.null(peaks)) {
      colnames(peaks) <- c(basenames, verbosenames)
      colnames(peakinfo) <- c("scale", "scaleNr", "scpos", "scmin", "scmax")
      for (p in 1:dim(peaks)[1]) {
        if (integrate == 1) {
          lm <- xcms:::descendMin(wCoefs[, peakinfo[p, "scaleNr"]], 
                                  istart = peakinfo[p, "scpos"])
          gap <- all(d[lm[1]:lm[2]] == 0)
          if ((lm[1] == lm[2]) || gap)
            lm <- xcms:::descendMinTol(d, startpos = c(peakinfo[p, "scmin"], 
                                                       peakinfo[p, "scmax"]), 
                                       maxDescOutlier)
        } else {
          lm <- xcms:::descendMinTol(d, startpos = c(peakinfo[p, "scmin"], 
                                                     peakinfo[p, "scmax"]), 
                                     maxDescOutlier)
        }
        lm <- xcms:::.narrow_rt_boundaries(lm, d)
        lm_seq <- lm[1]:lm[2]
        pd <- d[lm_seq]
        peakrange <- td[lm]
        peaks[p, "rtmin"] <- scantime[peakrange[1]]
        peaks[p, "rtmax"] <- scantime[peakrange[2]]
        peaks[p, "maxo"] <- max(pd)
        pwid <- (scantime[peakrange[2]] - 
                   scantime[peakrange[1]])/(peakrange[2] - peakrange[1])
        if (is.na(pwid)) 
          pwid <- 1
        peaks[p, "into"] <- pwid * sum(pd)
        db <- pd - baseline
        peaks[p, "intb"] <- pwid * sum(db[db > 0])
        peaks[p, "lmin"] <- lm[1]
        peaks[p, "lmax"] <- lm[2]
        if (fitgauss) {
          td_lm <- td[lm_seq]
          md <- max(pd)
          d1 <- pd/md
          pgauss <- xcms:::fitGauss(td_lm, pd, 
                                    pgauss = list(mu = peaks[p, "scpos"], 
                                                  sigma = peaks[p, "scmax"] - 
                                                    peaks[p, "scmin"], 
                                                  h = peaks[p, "maxo"]))
          rtime <- peaks[p, "scpos"]
          if (!any(is.na(pgauss)) && all(pgauss > 0)) {
            gtime <- td[match(round(pgauss$mu), td)]
            if (!is.na(gtime)) {
              rtime <- gtime
              peaks[p, "mu"] <- pgauss$mu
              peaks[p, "sigma"] <- pgauss$sigma
              peaks[p, "h"] <- pgauss$h
              peaks[p, "egauss"] <- sqrt((1/length(td_lm)) * 
                                           sum(((d1 - gauss(td_lm, pgauss$h/md, 
                                                            pgauss$mu, pgauss$sigma))^2)))
            }
          }
          peaks[p, "rt"] <- scantime[rtime]
          if (peaks[p, "rt"] < peaks[p, "rtmin"]) 
            peaks[p, "rt"] <- scantime[peaks[p, "scpos"]]
        }
        else peaks[p, "rt"] <- scantime[peaks[p, "scpos"]]
      }
      peaks <- xcms:::joinOverlappingPeaks(td, d, otd, omz, od, 
                                           scantime, sr, peaks, maxGaussOverlap = 0.5,
                                           mzCenterFun = mzCenterFun)
    }
    if (!is.null(peaks)) {
      peaklist[[length(peaklist) + 1]] <- peaks
    }
  }
  close(pb)
  p <- do.call(rbind, peaklist)
  if (!verboseColumns) 
    p <- p[, basenames, drop = FALSE]
  uorder <- order(p[, "into"], decreasing = TRUE)
  pm <- as.matrix(p[, c("mzmin", "mzmax", "rtmin", 
                        "rtmax"), drop = FALSE])
  uindex <- xcms:::rectUnique(pm, uorder, mzdiff, ydiff = -1e-05)
  pr <- p[uindex, , drop = FALSE]
  message(" OK: ", nrow(pr), " found.")
  return(pr)
}

.milliWave(mz, int, scantime, valsPerSpect, ppm, peakwidth, snthresh, prefilter,
           mzCenterFun, integrate, mzdiff, fitgauss, noise, sleep, 
           verboseColumns, roiList, firstBaselineCheck, roiScales)
