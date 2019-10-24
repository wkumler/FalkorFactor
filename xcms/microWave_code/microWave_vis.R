
# Setup things ----
library(dplyr)

source("xcms/microWave_code/microWave_functions.R")
load("xcms/microWave_code/temp_eic_list")
peak_df <- read.csv("xcms/microWave_code/temp_peak_df.csv")

peak_df <- mutate(peak_df, qscore=Peak_SNR*Peak_gauss_fit^4*log10(Peak_area))
# BUILD THIS NEXT STEP INTO PEAKFINDING STEP AND CHECK FOR BUGS
peak_df <- mutate(peak_df, Peak_id=as.character(Peak_id))
peak_df <- arrange(peak_df, desc(qscore))

peakCheck(eic_list, peak_df, "1.1.1")
# for(i in peak_df$Peak_id){
#   peakCheck(eic_list, peak_df, i)
#   replot <- readline(prompt = "Press Enter") 
#   if(replot=="j"){
#     peakCheck(eic_list, peak_df, i, zoom=T)
#   } else if(replot=="k") {
#     peakCheck(eic_list, peak_df, i, pts = T, zoom = T)
#   } else if(replot=="m"){
#     v <- getMetlinMz(peak_df[peak_df$Peak_id==i,"Peak_mz"]-1.007825)
#   } else {
#     next
#   }
#   readline(prompt = "Continue?")
# }

peak_df <- filter(peak_df, qscore>1)
peak_df_best <- filter(peak_df, qscore>5)


renderPeakOverview(peak_df_best)
#pdf(file = "xcms/microWave_code/TempPeakPlot.pdf")
#renderPeakOverview(peak_df_best)
#dev.off()


# Isotopes!
peak_df_best$Isotopes <- "None"
for(i in seq_len(nrow(peak_df_best))){
  given_isos <- findIsos(given_peak_id = peak_df_best[i,"Peak_id"], 
                         peak_df = peak_df, eic_list = eic_list)
  if(nrow(given_isos)>0){
    peak_df_best$Isotopes[i] <- list(split(given_isos[c("Peak_id", "mz_match", "rt_match", "cor")], seq_len(nrow(given_isos))))
  }
}
peak_df_best[peak_df_best$Isotopes!="None",c("Peak_id", "Isotopes")]

isoCheck(peak_df = peak_df, eic_list = eic_list, peak_id_1 = "1.1.22", "3.4.5")
