# Testing new scale generation algorithm

roi <- roi_list[[2]]@data
lmaoPlotEm(roi)
diagnoseWavelets(roi)

roi_intensity <- as.numeric(table(cut(rnorm(10000), breaks = 20)))
roi_intensity <- c(roi_intensity, rep(0, 1024-length(roi_intensity)))
plot(roi_intensity, type="l")
wCoefs <- xcms:::MSW.cwt(roi_intensity, scales)
plot(roi_intensity, type="l", lwd=2)
for(i in seq_along(colnames(wCoefs))){
  lines(wCoefs[,i], col=rainbow(length(colnames(wCoefs)))[i])
}

scales <- round(seq(1, 2^ceiling(log2(length(roi_intensity)))/12, length.out = 11))
scales <- seq(11,44)

possible_peaks <- xcms:::MSW.cwt(roi$int, scales, wavelet = "mexh") %>%
  xcms:::MSW.getLocalMaximumCWT() %>%
  xcms:::MSW.getRidge()








  
# Legible xcms:::MSW.cwt code below
# Extend the ROI to an FFT-acceptable length (64, 128, 256, etc.)
# Done by mirroring the last [however many scans] which are later removed
roi_intensity_extended <- xcms:::MSW.extendNBase(roi_intensity, nLevel = NULL, base = 2)


# xcms:::MSW.cwt
# Define wavelet shape
wavelet_xvals <- seq(-6, 6, length.out = 256)
wavelet_shape <- (2/sqrt(3) * pi^(-0.25)) * (1 - wavelet_xvals^2) * exp(-wavelet_xvals^2/2)




psi_xval <- seq(0,12, length.out = 256)
dxval <- psi_xval[2]
xmax <- 12

wCoefs <- list()
for (i in 1:length(scales)) {
  scale.i <- scales[i]
  mex_h <- rep(0, length(roi_intensity_extended))
  # WHAT DECIDES THE LENGTH OF THIS LINE
  sample_points <- 1 + floor((0:(scale.i * 12))/(scale.i * 12/255))
  if (length(sample_points) == 1) 
    sample_points <- c(1, 1)
  mex_h[1:length(sample_points)] <- rev(wavelet_shape[sample_points]) - mean(wavelet_shape[sample_points])
  if (length(mex_h) > length(roi_intensity_extended)) {
    i <- i - 1
    break
  }
  wCoefs.i <- 1/sqrt(scale.i) * convolve(roi_intensity_extended, mex_h)
  wCoefs.i <- c(wCoefs.i[(length(roi_intensity_extended) - floor(length(sample_points)/2) + 1):length(roi_intensity_extended)], 
                wCoefs.i[1:(length(roi_intensity_extended) - floor(length(sample_points)/2))])
  wCoefs[[i]] <- wCoefs.i
}
wCoefs <- do.call(cbind, wCoefs)
if (i < 1) # If there were no scales that magically worked?
  return(NA)
scales <- scales[1:i]
if (length(scales) == 1) 
  wCoefs <- matrix(wCoefs, ncol = 1)
colnames(wCoefs) <- scales
wCoefs <- wCoefs[1:length(roi_intensity), , drop = FALSE]

par(mfrow=c(4,1))
plot(roi_intensity_extended, type = "l")
plot(mex_h, type = "l")
#plot(convolve(roi_intensity_extended, mex_h), type = "l")
plot(convolve(roi_intensity_extended, mex_h[c((length(roi_intensity_extended) - floor(length(sample_points)/2) + 1):length(roi_intensity_extended),
                          1:(length(roi_intensity_extended) - floor(length(sample_points)/2)))]), 
     type = "l")

plot(wCoefs[,ncol(wCoefs)], type="n")
for(i in ncol(wCoefs):1){
  lines(wCoefs[,i], col=rainbow(ncol(wCoefs))[i], lwd=2)
}
ncol(wCoefs)
colnames(wCoefs)
