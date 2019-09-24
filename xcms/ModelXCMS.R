library(xcms)
source("ModelXCMS_source.R")


model_peak <- as.numeric(table(cut(rnorm(20000), breaks = 40)))

# Working with a model peak
long_scans <- 1:1000
long_ints <- c(numeric(480), model_peak, numeric(480))
long_chrom <- Chromatogram(rtime = long_scans, intensity = long_ints, mz = 250)
plot(long_chrom)
# Not found by default because mapply simplifies
findChromPeaks(long_chrom, param = CentWaveParam())
# Works if SIMPLIFY=F is added to mapply call within .getRtROI
newmapply_findChromPeaks(long_chrom, param = CentWaveParam())


# Working with a shorter peak to illustrate SNR problems
short_scans <- 201:400
short_ints <- c(numeric(80), model_peak, numeric(80))
short_chrom <- Chromatogram(rtime = short_scans, intensity = short_ints, mz = 250)
plot(short_chrom)
# Not found because the SNR is messed up
newmapply_findChromPeaks(short_chrom, param = CentWaveParam())
# SNR must be disabled to see peak
newmapply_findChromPeaks(short_chrom, param = CentWaveParam(snthresh = 1))
# Alternatively, adding a tiny amount of noise also resolves the problem
short_chrom@intensity <- short_chrom@intensity + runif(200)
newmapply_findChromPeaks(short_chrom, param = CentWaveParam())


# Working with an even shorter peak to illustrate ROI length <-> wavelets dependence
shortest_scans <- 201:250
shortest_ints <- c(numeric(5), model_peak, numeric(5)) + round(runif(50), digits = 2)
shortest_chrom <- Chromatogram(rtime = shortest_scans, intensity = shortest_ints, mz = 250)
plot(shortest_chrom)
newmapply_findChromPeaks(shortest_chrom, param = CentWaveParam(snthresh = 1), fixMSW=F)
newmapply_findChromPeaks(shortest_chrom, param = CentWaveParam(snthresh = 1), fixMSW=T)

