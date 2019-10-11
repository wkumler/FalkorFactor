
# Which wavelet is most similar in shape to the peak?

peak <- as.numeric(table(cut(rnorm(100000), breaks = 100)))
peak_ext <- c(peak, integer(2^(ceiling(log2(50*12)))-length(peak)))
wcoefs_ext <- xcms:::MSW.cwt(peak_ext, scales = 1:50)
wcoefs <- wcoefs_ext[1:length(peak),]

plot(peak, type = "l", lwd=2, ylim=c(min(wcoefs), max(wcoefs)))
apply(wcoefs, 2, lines)

apply(wcoefs, 2, cor, y=peak)



layout(matrix(1:49, ncol=7, byrow = T))
par(mar=c(0,0,0,0))
for(i in seq_len(ncol(wcoefs))){
  plot(wcoefs[,i], peak, axes=F)
  model_i <- lm(peak~wcoefs[,i])
  abline(model_i)
  legend("topleft", legend = round(summary(model_i)$r.squared, 2))
}
layout(1)




widthFinder <- function(peak_ints, peak_center, peakwidth){
  peak_widths_to_check <- seq(min(peakwidth), max(peakwidth), 2)
  peak_ints_buffered <- c(numeric(max(peakwidth)/2), peak_ints, numeric(max(peakwidth)/2))
  peak_fits <- sapply(peak_widths_to_check, function(pred_peak_width){
    perf_peak <- exp((-seq(-2, 2, length.out = pred_peak_width+1)^2))
    peak_left <- (peak_center+max(peakwidth)/2-pred_peak_width/2)
    peak_right <- (peak_center+max(peakwidth)/2+pred_peak_width/2)
    relevant_ints <- peak_ints_buffered[peak_left:peak_right]
    return(cor(relevant_ints, perf_peak))
  })
  return(c(peak_widths_to_check[which.max(peak_fits)]))
}



peak <- as.numeric(table(cut(rnorm(100000), breaks = 50)))
peak <- as.numeric(table(cut(rbeta(10000, 100, 4), breaks = 80)))
plot(peak/max(peak), type="b")
peak_width <- widthFinder(peak_ints = peak, 
                         peak_center = which.max(peak), 
                         peakwidth = c(20, 80))
print(peak_width)
plot(peak/max(peak), type="b")
perf_peak <- exp((-seq(-2, 2, length.out = peak_width[1]+1)^2))
id_idx <- floor(which.max(peak)-peak_width[1]/2):ceiling(which.max(peak)+peak_width[1]/2)
polygon(x = c(id_idx[1], id_idx, id_idx[length(id_idx)]), 
        y = c(0, peak[id_idx]/max(peak), 0), col = "#FF000022")
