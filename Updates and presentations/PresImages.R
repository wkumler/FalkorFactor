load("xcms/raw_data")
x <- raw_data %>%
  filterMsLevel(msLevel. = 1L) %>%
  selectFeatureData(fcol = c(MSnbase:::.MSnExpReqFvarLabels, "centroided")) %>%
  lapply(1:length(fileNames(.)), FUN=filterFile, object = .) %>%
  `[[`(1) %>%
  spectra()

mzs <- lapply(x, mz)
mz <- unlist(mzs, use.names = FALSE)
int <- unlist(lapply(x, intensity), use.names = FALSE)
rt <- unlist(lapply(x, rtime), use.names = FALSE)
rts <- rep(rt, sapply(mzs, length))

all_data <- data.frame(mz, int, rts)
scales <- seq(1, 40, by=4)

roi <- all_data %>% 
  filter(mz>90.0910&mz<90.0925)

wcoef_matrix <- xcms:::MSW.cwt(roi$int, scales = scales, wavelet = "mexh")
local_maxima <- xcms:::MSW.getLocalMaximumCWT(wcoef_matrix)

roi <- roi_list[[4]]

roi <- mutate(roi, mrt=rt/60)
ploy <- filter(roi, mrt>8&mrt<8.75)

png(filename = "Slideback.png", width = 1920, height = 1080, units = "px", type = "cairo")
plot(roi$mrt, roi$int, type="l", lwd=6, xlim=c(7.95, 9), axes=F, ylab="", xlab="")
lines(ploy$mrt, ploy$int, col="red", lwd=6)
dev.off()

png(filename = "Samplepeak.png", width = 1920, height = 1080, units = "px", type = "cairo")
plot(roi$mrt, roi$int, type="l", lwd=6, xlim=c(7.95, 9), axes=F, ylab="", xlab="")
x <- c(ploy$mrt[1], ploy$mrt, ploy$mrt[length(ploy$mrt)])
y <- c(0, ploy$int, 0)
polygon(x, y, col = rgb(1,0,0,0.2), border = "red", lwd=6)
dev.off()

png(filename = "FittedPeak.png", width = 1920, height = 1080, units = "px", type = "cairo")
plot(roi$rt, roi$int, type="l", xlim=c(100, 800), lwd=6, xaxt="n", yaxt="n", ylab="", xlab="")
yvals <- ((wcoef_matrix[,9]+80000)/max(wcoef_matrix[,9]))*max(roi$int)
lines(roi$rt, yvals, col=rainbow(10)[4], lwd=10)
dev.off()


scene = list(camera = list(eye = list(x = -2, y = 0.1, z = 0.5)),
             yaxis = list(autorange = "reversed", title = "Retention time"),
             xaxis = list(title = "Wavelet scale", 
                          ticktext = list(colnames(wcoef_matrix)), 
                          tickvals = list(seq(0, 9, length.out = ncol(wcoef_matrix))),
                          tickmode = "array"),
             zaxis = list(title = "Wavelet intensity"))
p <- plotly::plot_ly(z = ~wcoef_matrix[225:850,], hoverinfo='skip') %>%
  plotly::add_surface(hoverinfo='skip', showscale=F, showlegend=F) %>%
  plotly::layout(title = "3D rendering of wavelet coefficients", 
                 scene = scene)
plotly::api_create(p, filename = "hover-text")
Sys.setenv("plotly_username"="wkumler")
Sys.setenv("plotly_api_key"="oijnX1fd4HWB9Uv9e7Up")

png(filename = "CleanPeak.png", width = 1920, height = 1080, units = "px", type = "cairo")
peak <- peak_df[3,]
roi <- roi_list[[peak$ROI_number]]
plot(roi$rt, roi$int, type="l", lwd=6, axes=F, ylab="", xlab="")
dev.off()

png(filename = "MessyPeak.png", width = 1920, height = 1080, units = "px", type = "cairo")
peak <- peak_df[234,]
roi <- roi_list[[peak$ROI_number]]
plot(roi$rt, roi$int, type="l", lwd=6, axes=F, ylab="", xlab="")
peak_roi <- filter(roi, rt>=peak$peak_lefts&rt<=peak$peak_rights)
lines(peak_roi$rt, peak_roi$int, col="red", lwd=6)
dev.off()

png(filename = "MessyPeak2.png", width = 1920, height = 270, units = "px", type = "cairo")
peak <- peak_df[236,]
roi <- roi_list[[peak$ROI_number]]
plot(roi$rt, roi$int, type="l", lwd=6, axes=F, ylab="", xlab="")
peak_roi <- filter(roi, rt>=peak$peak_lefts&rt<=peak$peak_rights)
lines(peak_roi$rt, peak_roi$int, col="red", lwd=6)
dev.off()


par(mfrow=c(2,10))
par(mar=c(0.1, 0.1, 0.1, 0.1))

coef_peak_df <- as.data.frame(do.call(rbind, peak_list)) %>% 
  mutate(qty=coef_areas) %>%
  arrange(desc(qty))
png(filename = "PeakCoefs.png", width = 960, height = 360, units = "px", type = "cairo")
par(mfrow=c(2,10))
par(mar=c(0.1, 0.1, 0.1, 0.1))
for(i in 1:20){
  peak <- coef_peak_df[i,]
  roi <- roi_list[[peak$ROI_number]]
  plot(roi$rt, roi$int, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
  peak_roi <- filter(roi, rt>=peak$peak_lefts&rt<=peak$peak_rights)
  lines(peak_roi$rt, peak_roi$int, col="red")
}
dev.off()

sharpness_peak_df <- as.data.frame(do.call(rbind, peak_list)) %>% 
  mutate(qty=ROI_sharpness) %>%
  arrange(desc(qty))
png(filename = "PeakSharpness.png", width = 960, height = 270, units = "px", type = "cairo")
par(mfrow=c(2,10))
par(mar=c(0.1, 0.1, 0.1, 0.1))
for(i in 1:20){
  peak <- sharpness_peak_df[i,]
  roi <- roi_list[[peak$ROI_number]]
  plot(roi$rt, roi$int, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
  peak_roi <- filter(roi, rt>=peak$peak_lefts&rt<=peak$peak_rights)
  lines(peak_roi$rt, peak_roi$int, col="red")
}
dev.off()

peak_height_df <- as.data.frame(do.call(rbind, peak_list)) %>% 
  mutate(qty=peak_tops) %>%
  arrange(desc(qty))
png(filename = "PeakHeight.png", width = 960, height = 270, units = "px", type = "cairo")
par(mfrow=c(2,10))
par(mar=c(0.1, 0.1, 0.1, 0.1))
for(i in 1:20){
  peak <- peak_height_df[i,]
  roi <- roi_list[[peak$ROI_number]]
  plot(roi$rt, roi$int, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
  peak_roi <- filter(roi, rt>=peak$peak_lefts&rt<=peak$peak_rights)
  lines(peak_roi$rt, peak_roi$int, col="red")
}
dev.off()

best_peak_df <- as.data.frame(do.call(rbind, peak_list)) %>% 
  mutate(qty=ROI_sharpness*coef_areas*peak_tops) %>%
  arrange(desc(qty))
png(filename = "PeakBest.png", width = 960, height = 270, units = "px", type = "cairo")
par(mfrow=c(2,10))
par(mar=c(0.1, 0.1, 0.1, 0.1))
for(i in 1:20){
  peak <- best_peak_df[i,]
  roi <- roi_list[[peak$ROI_number]]
  plot(roi$rt, roi$int, type="l", xaxt="n", yaxt="n", xlab="", ylab="")
  peak_roi <- filter(roi, rt>=peak$peak_lefts&rt<=peak$peak_rights)
  lines(peak_roi$rt, peak_roi$int, col="red")
}
dev.off()



par(mfrow=c(2,2))
par(mar=c(1.1, 1.1,1.1,1.1))
plot(sort(cumsum(rnorm(100, mean = 1))), type="l", xaxt="n", yaxt="n", ylab="", xlab="", lwd=2)
plot(rev(sort(cumsum(rnorm(100, mean = 1)))), type="l", xaxt="n", yaxt="n", ylab="", xlab="", lwd=2)
plot((1:100)^2, type="l", xaxt="n", yaxt="n", ylab="", xlab="", lwd=2)
plot(1:100, type="l", xaxt="n", yaxt="n", ylab="", xlab="", lwd=2)