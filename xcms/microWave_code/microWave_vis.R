
# Setup things ----
library(dplyr)

source("xcms/microWave_code/microWave_functions.R")


peak_df <- arrange(peak_df, desc(qscore))

#Look at all the pretty peaks! ----
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



# Visualize the isotopes
isotope_df <- filter(peak_df, !is.na(Isotopes)) #%>% select("Peak_id", "Isotopes")
isoCheck(peak_df = peak_df, eic_list = eic_list, 
         peak_id_1 = isotope_df$Peak_id[1], 
         peak_id_2 = isotope_df$Isotopes[[1]]$`1`$Peak_id)
# for(i in seq_len(nrow(isotope_df))){
#   isoCheck(peak_df = peak_df, eic_list = eic_list,
#            peak_id_1 = isotope_df$Peak_id[i],
#            peak_id_2 = isotope_df$Isotopes[[i]]$`1`$Peak_id)
# }



# Decide which ones are worth keeping (quality score>1)
# and which ones are worth analyzing (quality score > 5)
peak_df <- filter(peak_df, qscore>1)
peak_df_best <- filter(peak_df, qscore>5)

# Visualize the highest-quality peaks
renderPeakOverview(peak_df_best)
# # Option to export to PDF
# pdf(file = "xcms/microWave_code/TempPeakPlot.pdf")
# renderPeakOverview(peak_df_best)
# dev.off()
