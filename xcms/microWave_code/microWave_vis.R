
# Setup things ----
library(dplyr)

source("xcms/microWave_code/microWave_functions.R")

peak_df <- arrange(peak_df, desc(qscore))

#Look at all the pretty peaks! ----
peakCheck(peak_df$Peak_id[1])
# for(i in peak_df$Peak_id){
#   peakCheck(i)
#   replot <- readline(prompt = "Press Enter")
#   if(replot=="j"){
#     peakCheck(i, zoom=T)
#   } else if(replot=="k") {
#     peakCheck(i, pts = T, zoom = T)
#   } else if(replot=="m"){
#     v <- getMetlinMz(peak_df[peak_df$Peak_id==i,"Peak_mz"]-1.007825)
#   } else {
#     next
#   }
#   readline(prompt = "Continue?")
# }



# Visualize the isotopes
isotope_peak_df <- filter(peak_df, !is.na(Isotopes)) %>% select("Peak_id", "Isotopes")
isoCheck(peak_id_1 = isotope_peak_df$Peak_id[1],
         peak_id_2 = isotope_peak_df$Isotopes[[1]]$`1`$Peak_id)
for(i in seq_len(nrow(isotope_peak_df))){
  isoCheck(peak_id_1 = isotope_peak_df$Peak_id[i],
           peak_id_2 = isotope_peak_df$Isotopes[[i]]$`1`$Peak_id, default_layout = F)
  #readline(prompt = "Press Enter")
}

# Check weird isotopes against isotopically-labeled standards
stds <- read.csv(paste0("https://raw.githubusercontent.com/kheal/Example_Unta",
                        "rgeted_Metabolomics_Workflow/master/Ingalls_Lab_Stan",
                        "dards.csv"), stringsAsFactors = F)
iso_stds <- subset(stds, grepl("D|C\\(13\\)", stds$Emperical.Formula)&stds$Column=="HILIC")
iso_stds_df <- data.frame(name=iso_stds$Compound.Name, 
                          formula=iso_stds$Emperical.Formula,
                          mzs=as.numeric(iso_stds$m.z),
                          rts=iso_stds$RT..min.*60)
for(i in seq_len(nrow(iso_stds_df))){
  std_mz <- iso_stds_df$mzs[i]
  std_info <- filter(peak_df, Peak_mz<1.000025*std_mz&Peak_mz>0.999975*std_mz)
  if(nrow(std_info)>0){
    print(cbind(iso_stds_df[i,], select(std_info, c(Peak_id, qscore, Isotopes))))
  }
}


# Decide which ones are worth keeping (quality score>5)
peak_df_best <- filter(peak_df, qscore>2)

# Visualize the highest-quality peaks
renderPeakOverview(peak_df_best)
# Option to export to PDF
# pdf(file = "xcms/microWave_code/TempPeakPlot.pdf")
# renderPeakOverview(peak_df_best)
# dev.off()
