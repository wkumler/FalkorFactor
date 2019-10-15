
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


possible_peakwidths <- min(peakwidth_scans):max(peakwidth_scans)
perf_peak_list <- lapply(possible_peakwidths, function(x){
  exp((-seq(-2, 2, length.out = x+1)^2))
})
names(perf_peak_list) <- as.character(possible_peakwidths)


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
  
  return(list(edges=c(floor(peak_center-best_peak_width/2), 
                      ceiling(peak_center+best_peak_width/2)),
              cor=max(peak_fits),
              residuals=residuals))
}



peak <- as.numeric(table(cut(rnorm(1000), breaks = 80)))
peak_width <- widthFinder(peak_ints = peak, 
                         peak_center = which.max(peak), 
                         peakwidth = c(21, 88))
layout(matrix(c(1,1,2), nrow=3))
plot(peak, type="b")
peak_background <- mean(c(head(peak, 3), tail(peak, 3)))
peak_noise <- sd(peak_width$residuals)*max(peak)
abline(h=peak_background+peak_noise)
legend("topleft", legend=peak_width$edges[2]-peak_width$edges[1])
legend("topright", legend=max(peak)/(peak_background+peak_noise))
plot(peak_width$residuals)