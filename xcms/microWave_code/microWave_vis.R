
library(dplyr)

source("xcms/microWave_code/microWave_functions.R")
load("xcms/microWave_code/temp_eic_list")
peak_df <- read.csv("xcms/microWave_code/temp_peak_df.csv")

peak_df <- mutate(peak_df, qscore=Peak_SNR*Peak_gauss_fit^4*log10(Peak_area))
peak_df <- arrange(peak_df, desc(qscore))

peakCheck(eic_list, peak_df, "1.1.1")
for(i in peak_df$Peak_id){
  peakCheck(eic_list, peak_df, i)
  replot <- readline(prompt = "Press Enter") 
  if(replot=="j"){
    peakCheck(eic_list, peak_df, i, zoom=T)
  } else if(replot=="k") {
    peakCheck(eic_list, peak_df, i, pts = T, zoom = T)
  } else if(replot=="m"){
    v <- getMetlinMz(peak_df[peak_df$Peak_id==i,"Peak_mz"]-1.007825)
  } else {
    next
  }
  readline(prompt = "Continue?")
}

peak_df_best <- filter(peak_df, qscore>3)

pdf(file = "xcms/microWave_code/TempPeakPlot.pdf")
ylimits <- c(min(peak_df_best$Peak_mz), max(peak_df_best$Peak_mz))
xlimits <- c(min(peak_df_best$Peak_start_time), max(peak_df_best$Peak_end_time))
plot(1, ylim=ylimits, xlim=xlimits)
factored_peakareas <- factor(round(log2(peak_df_best$Peak_area)))
peak_shades <- gray.colors(length(unique(factored_peakareas)), start = 0, end = 1)[factored_peakareas]
peak_shades <- hcl.colors(length(unique(factored_peakareas)), rev = T)[factored_peakareas]
for(i in seq_len(nrow(peak_df_best))){
  segments(x0=peak_df_best[i, "Peak_start_time"], x1=peak_df_best[i, "Peak_end_time"],
           y0=peak_df_best[i, "Peak_mz"], y1=peak_df_best[i, "Peak_mz"],
           col = peak_shades[i], lwd=floor(log(peak_df_best[i, "qscore"])))
  text(x = mean(c(peak_df_best[i, "Peak_start_time"], peak_df_best[i, "Peak_end_time"])),
       y = peak_df_best[i, "Peak_mz"]+0.2, labels = peak_df_best[i, "Peak_id"], 
       cex = 0.5, col = peak_shades[i])
}
abline(h = floor(min(peak_df$Peak_mz)):ceiling(max(peak_df$Peak_mz)), lty=2, col="gray")
dev.off()

# Isotopes between 100 and 120 Da
isoCheck(peak_df, eic_list, "1.1.1", "3.4.1")
isoCheck(peak_df, eic_list, "12.1.4", "116.4.1")
isoCheck(peak_df, eic_list, "5.1.2", "57.1.1")
isoCheck(peak_df, eic_list, "43.4.1", "60.1.1")
isoCheck(peak_df, eic_list, "7.1.7", "81.1.1")
isoCheck(peak_df, eic_list, "7.1.2", "81.2.1")


isoCheck(peak_df, eic_list, "5.1.1", "40.2.1")
isoCheck(peak_df, eic_list, "37.1.2", "15.1.2")
isoCheck(peak_df, eic_list, "17.1.1", "115.1.1")
