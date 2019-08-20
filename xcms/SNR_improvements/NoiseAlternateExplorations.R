# xcms:::getLocalNoiseEstimate
load("xcms/SNR_improvements/NoiseAlternateWorkspace")

Nscantime = length(rt)
threshold=noise
num=minPtsAboveBaseLine

if (length(d) < Nscantime) {
  drange <- which(td %in% ftd)
  n1 <- d[-drange]
  n1.cp <- continuousPtsAboveThresholdIdx(n1, threshold = threshold, 
                                          num = num)
  n1 <- n1[!n1.cp]
  if (length(n1) > 1) {
    baseline1 <- mean(n1)
    sdnoise1 <- sd(n1)
  }
  else baseline1 <- sdnoise1 <- 1
  d1 <- drange[1]
  d2 <- drange[length(drange)]
  nrange2 <- c(max(1, d1 - noiserange[1]):d1, d2:min(length(d), 
                                                     d2 + noiserange[1]))
  n2 <- d[nrange2]
  n2.cp <- continuousPtsAboveThresholdIdx(n2, threshold = threshold, 
                                          num = num)
  n2 <- n2[!n2.cp]
  if (length(n2) > 1) {
    baseline2 <- mean(n2)
    sdnoise2 <- sd(n2)
  }
  else baseline2 <- sdnoise2 <- 1
} else {
  trimmed <- xcms:::trimm(d, c(0.05, 0.95))
  baseline1 <- baseline2 <- mean(trimmed)
  sdnoise1 <- sdnoise2 <- sd(trimmed)
}
c(min(baseline1, baseline2), min(sdnoise1, sdnoise2))


x <- d
trim = c(0.05, 0.95)
a <- sort(x[x > 0])
Na <- length(a)
quant <- round((length(a) * trim[1]) + 1):round(length(a) * trim[2])
xcms_noise_vals <- a[quant]
will_noise_vals <- head(a, -max(peakwidth))
sd_noise_vals <- a[a<(mean(a)+sd(a))]
iqr_noise_vals <- a[a<(median(a)+IQR(a))]




x <- as.data.frame(d) %>%
  mutate(index=1:length(d)) %>%
  mutate(xcms_noise_vals=d%in%xcms_noise_vals[round((length(d) * trim[1]) + 1):round(length(d) * trim[2])]) %>%
  mutate(will_noise_vals=d%in%head(sort(d), -max(peakwidth))) %>%
  mutate(sd_noise_vals=d<(mean(d)+sd(d))) %>%
  mutate(iqr_noise_vals=d<(median(d)+IQR(d)))
toPlot <- function(what, color, main) {
  plot(x$index, x$d, main=main)
  w <- filter(x, x[[what]])
  points(w$index, w$d, col=color)
  legend("topleft", legend = paste("SD:", round(sd(w$d))))
  boxplot(w$d, width = 1, axes=F, cex=1)
}


layout(matrix(c(1,1,1,2,3,3,3,4,5,5,5,6,7,7,7,8,9,9,9,10), nrow = 5, byrow = T))
par(mar=c(2.1, 0.1, 2.1, 0.1))

plot(a, main="Initial peak values")
abline(v=0, col="blue")
abline(v=0.5, col="green")
abline(v=length(sd_noise_vals), col="blue")
abline(v=which(a==min(xcms_noise_vals)), col="red")
abline(v=which(a==max(xcms_noise_vals)), col="red")
abline(v=length(will_noise_vals), col="green")
abline(v=1, col="purple")
abline(v=length(iqr_noise_vals), col="purple")
legend("topleft", legend = paste("SD:", round(sd(a))))
boxplot(a, axes=F, width=1, cex=1)

toPlot("xcms_noise_vals", "red", "Top and bottom 5% removed")
toPlot("sd_noise_vals", "blue", "1 SD above mean removed")
toPlot("will_noise_vals", "green", "Top 'max peakwidth' removed")
toPlot("iqr_noise_vals", "purple", "Beyond IQR removed")
par(mfrow=c(1,1))








