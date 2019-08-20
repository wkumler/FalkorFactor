# Detailed review of all the things that happen in findPeaksCentWave



will_peaksWithCentWave <- function(int, rt, param){

#Setup things ----

load("xcms/chr_raw")


peakwidth = c(20, 50)
snthresh = 10
prefilter = c(3, 100)
integrate = 1
fitgauss = FALSE
noise = 0
verboseColumns = FALSE
firstBaselineCheck = TRUE

maxGaussOverlap <- 0.5

peaklist <- list()
peakinfo_names <- c("scale", "scaleNr", "scpos", "scmin", "scmax")
basenames <- c("mz", "mzmin", "mzmax", "rt", "rtmin", "rtmax", 
               "into", "intb", "maxo", "sn")
verbosenames <- c("egauss", "mu", "sigma", "h", "f", "dppm", 
                  "scale", "scpos", "scmin", "scmax", "lmin", "lmax")
peaks_names <- c(basenames, verbosenames)
peaks_ncols <- length(peaks_names)



# Startup checks ----
if (length(peakwidth) != 2) 
  stop("'peakwidth' has to be a numeric of length 2")
if(any(peakwidth<=0))
  stop("'peakwidth' must be a positive integer")
int[is.na(int)] <- 0



# Identify regions of interest by finding local maxima ----
plot(int, ylim=c(0, max(int)+1000))
rois <- xcms:::.getRtROI(int, rt, peakwidth = peakwidth, noise = noise, 
                         prefilter = prefilter)
for(row in 1:dim(rois)[1]){
  abline(v=rois[row, "sccent"], col="red")
  arrows(x0 = rois[row, "scmin"], x1 = rois[row, "scmax"],
         y0 = int[rois[row, "sccent"]], y1 = int[rois[row, "sccent"]],
         length = 0, col = "blue")
}



# Define scales ----
scalerange <- round((peakwidth/mean(diff(rt)))/2) #? Normalize peakwidth to distance between scans, then divide by 2?
scales <- seq(scalerange[1], scalerange[2])

minPeakWidth <- scales[1] # Min peak width is now not actually the minimum peak width, but instead about half?

noiserange <- c(min(scales) * 3, max(scales) * 3) #Will be used to walk a few peak widths outside to get a sense of the local noise
scRangeTol <- maxDescOutlier <- floor(minPeakWidth/2)
minCentroids <- minPtsAboveBaseLine <- max(4, minPeakWidth - 2)



# Run peak detection on each region of interest ----
for (i in seq_len(nrow(rois))) {
# i <- 3
  scmin <- rois[i, "scmin"]
  scmax <- rois[i, "scmax"]
  peaks <- matrix(ncol = peaks_ncols, nrow = 0, dimnames = list(character(), 
                                                                peaks_names))
  peakinfo <- matrix(ncol = 5, nrow = 0, dimnames = list(character(), 
                                                         peakinfo_names))
  
  td <- max(1, scmin - max(noiserange)):min(length(rt), scmax + max(noiserange))
  d <- int[td]
  scan.range <- c(min(td), max(td))
  otd <- scmin:scmax
  od <- int[otd]
  ftd <- max(td[1], scmin - scRangeTol):min(td[length(td)], scmax + scRangeTol)
  fd <- int[ftd]
  
  plot(d~td, main=paste("ROI number", i))
  points(fd~ftd, col = "blue")
  points(od~otd, col="red")
  
  
  
  # Estimate local noise
  if ((scmax - scmin + 1) >= 10 * minPeakWidth) {
    noised <- int
  } else {
    noised <- d
  }
  noise <- xcms:::estimateChromNoise(noised, trim = 0.05, 
                                     minPts = 3 * minPeakWidth)
  if (firstBaselineCheck & !xcms:::continuousPtsAboveThreshold(fd, threshold = noise, num = minPtsAboveBaseLine)){
    #print(paste("ROI number", i, "failed baseline check"))
    next #If we're doing a baseline check and this ROI fails it, skip to next ROI
  }
  
  source("xcms/NoiseAlternateFunctions.R")
  lnoise <- xcms:::getLocalNoiseEstimate(d, td, ftd, noiserange, 
                                         length(rt), threshold = noise, 
                                         num = minPtsAboveBaseLine)
  #lnoise <- getLocalNoise_maxpeaks(d, peakwidth = max(peakwidth))
  #lnoise <- getLocalNoise_sd(d)
  lnoise <- getLocalNoise_IQR(d)
  baseline <- max(1, min(lnoise[1], noise))
  sdnoise <- max(1, lnoise[2])
  sdthr <- sdnoise * snthresh
  if (!(any(fd - baseline >= sdthr))) {
    #print(paste("ROI number", i, "failed SNR check"))
    next
  }
  
  # Generate wavelets for remaining ROIs
  wCoefs <- xcms:::MSW.cwt(d, scales = scales, wavelet = "mexh")
  
  par(mfrow=c(2,1))
  par(mar=c(0.1, 0.1, 0.1, 0.1))
  plot(d, type="l", lwd=2, ylim=c(0, max(d)), axes=F)
  plot(wCoefs[,1], ylim=c(min(wCoefs), max(wCoefs)), axes=F)
  for(i in 1:dim(wCoefs)[2]){
    points(wCoefs[,i], col=rainbow(dim(wCoefs)[2])[i])
  }
  legend("topleft", legend = scales[1:dim(wCoefs)[2]], col = rainbow(dim(wCoefs)[2]), 
         lwd = 2, title = "Wavelet scale")
  par(mfrow=c(1,1))
  par(mar=c(4.1, 4.1, 0.1, 0.1))
  
  
  if (!(!is.null(dim(wCoefs)) && any((wCoefs - baseline) >= sdthr))) 
    next
  if (td[length(td)] == length(rt))
    wCoefs[nrow(wCoefs), ] <- wCoefs[nrow(wCoefs) - 1, 
                                     ] * 0.99
  localMax <- xcms:::MSW.getLocalMaximumCWT(wCoefs)
  rL <- xcms:::MSW.getRidge(localMax)
  wpeaks <- sapply(rL, function(x) {
    w <- min(1:length(x), ncol(wCoefs))
    any((wCoefs[x, w] - baseline) >= sdthr)
  })
  if (any(wpeaks)) {
    wpeaksidx <- which(wpeaks)
    for (p in 1:length(wpeaksidx)) {
      opp <- rL[[wpeaksidx[p]]]
      pp <- unique(opp)
      if (length(pp) >= 1) {
        dv <- td[pp] %in% ftd
        if (any(dv)) {
          if (any(d[pp[dv]] - baseline >= sdthr)) {
            inti <- numeric(length(opp))
            irange <- rep(ceiling(scales[1]/2), length(opp))
            for (k in 1:length(opp)) {
              kpos <- opp[k]
              r1 <- ifelse(kpos - irange[k] > 1, kpos - 
                             irange[k], 1)
              r2 <- ifelse(kpos + irange[k] < length(d), 
                           kpos + irange[k], length(d))
              inti[k] <- sum(d[r1:r2])
            }
            maxpi <- which.max(inti)
            if (length(maxpi) > 1) {
              m <- wCoefs[opp[maxpi], maxpi]
              bestcol <- which(m == max(m), arr.ind = TRUE)[2]
              best.scale.nr <- maxpi[bestcol]
            }
            else {
              best.scale.nr <- maxpi
            }
            best.scale <- scales[best.scale.nr]
            best.scale.pos <- opp[best.scale.nr]
            pprange <- min(pp):max(pp)
            lwpos <- max(1, best.scale.pos - best.scale)
            rwpos <- min(best.scale.pos + best.scale, 
                         length(td))
            p1 <- match(td[lwpos], otd)[1]
            p2 <- match(td[rwpos], otd)
            p2 <- p2[length(p2)]
            if (is.na(p1)) 
              p1 <- 1
            if (is.na(p2)) 
              p2 <- scmax - scmin + 1
            maxint <- max(od[p1:p2])
            peaks <- rbind(peaks, c(1, 1, 1, NA, NA, 
                                    NA, NA, NA, maxint, round((maxint - baseline)/sdnoise), 
                                    NA, NA, NA, NA, i, NA, best.scale, td[best.scale.pos], 
                                    td[lwpos], td[rwpos], NA, NA))
            peakinfo <- rbind(peakinfo, c(best.scale, 
                                          best.scale.nr, best.scale.pos, lwpos, 
                                          rwpos))
          }
        }
      }
    }
  }
  for (p in seq_len(nrow(peaks))) {
    if (integrate == 1) {
      lm <- xcms:::descendMin(wCoefs[, peakinfo[p, "scaleNr"]], 
                              istart = peakinfo[p, "scpos"])
      gap <- all(d[lm[1]:lm[2]] == 0)
      if ((lm[1] == lm[2]) || gap) 
        lm <- xcms:::descendMinTol(d, startpos = c(peakinfo[p, "scmin"], 
                                                   peakinfo[p, "scmax"]), 
                                   maxDescOutlier)
    }
    else {
      lm <- xcms:::descendMinTol(d, startpos = c(peakinfo[p, "scmin"], 
                                                 peakinfo[p, "scmax"]), 
                                 maxDescOutlier)
    }
    lm <- xcms:::.narrow_rt_boundaries(lm, d)
    lm_range <- lm[1]:lm[2]
    pd <- d[lm_range]
    peakrange <- td[lm]
    peaks[p, "rtmin"] <- rt[peakrange[1]]
    peaks[p, "rtmax"] <- rt[peakrange[2]]
    peaks[p, "maxo"] <- max(pd)
    pwid <- (rt[peakrange[2]] - rt[peakrange[1]])/(peakrange[2] - 
                                                     peakrange[1])
    if (is.na(pwid)) 
      pwid <- 1
    peaks[p, "into"] <- pwid * sum(pd)
    db <- pd - baseline
    peaks[p, "intb"] <- pwid * sum(db[db > 0])
    peaks[p, "lmin"] <- lm[1]
    peaks[p, "lmax"] <- lm[2]
    if (fitgauss) {
      td_lm <- td[lm_range]
      md <- max(pd)
      d1 <- pd/md
      pgauss <- fitGauss(td_lm, pd, pgauss = list(mu = peaks[p, 
                                                             "scpos"], sigma = peaks[p, "scmax"] - peaks[p, 
                                                                                                         "scmin"], h = peaks[p, "maxo"]))
      rtime <- peaks[p, "scpos"]
      if (!any(is.na(pgauss)) && all(pgauss > 0)) {
        gtime <- td[match(round(pgauss$mu), td)]
        if (!is.na(gtime)) {
          rtime <- gtime
          peaks[p, "mu"] <- pgauss$mu
          peaks[p, "sigma"] <- pgauss$sigma
          peaks[p, "h"] <- pgauss$h
          peaks[p, "egauss"] <- sqrt((1/length(td_lm)) * 
                                       sum(((d1 - gauss(td_lm, pgauss$h/md, pgauss$mu, 
                                                        pgauss$sigma))^2)))
        }
      }
      peaks[p, "rt"] <- rt[rtime]
      if (peaks[p, "rt"] < peaks[p, "rtmin"]) 
        peaks[p, "rt"] <- rt[peaks[p, "scpos"]]
    }
    else {
      peaks[p, "rt"] <- rt[peaks[p, "scpos"]]
    }
  }
  peaks <- xcms:::joinOverlappingPeaks(td, d, otd, rep(1, length(otd)), 
                                       od, rt, scan.range, peaks, maxGaussOverlap, mzCenterFun = mzCenter.wMean)
  if (!is.null(peaks)) 
    peaklist[[length(peaklist) + 1]] <- peaks
} # Goes with for loop at top
if (length(peaklist) == 0) {
  warning("No peaks found!")
  if (verboseColumns) 
    nopeaks <- matrix(nrow = 0, ncol = peaks_ncols, dimnames = list(character(), 
                                                                    peaks_names))
  else nopeaks <- matrix(nrow = 0, ncol = length(basenames), 
                         dimnames = list(character(), basenames))
  return(nopeaks[, -(1:3), drop = FALSE])
}
p <- do.call(rbind, peaklist)
if (!verboseColumns) 
  p <- p[, basenames, drop = FALSE]
uorder <- order(p[, "into"], decreasing = TRUE)
pm <- p[, c("mzmin", "mzmax", "rtmin", "rtmax"), drop = FALSE]
uindex <- xcms:::rectUnique(pm, uorder, xdiff = 0, ydiff = -1e-05)
unique(p[uindex, -(1:3), drop = FALSE])
}
