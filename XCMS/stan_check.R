library(tidyverse)
library(plotly)
library(data.table)

pmppm <- function(mass, ppm=4){
  if(mass<200){
    as.numeric(mass)+(c(-ppm, ppm)*200/1000000)
  } else {
    c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))
  }
}


polarity <- "pos"
pretty_folder <- paste0(polarity, "_pretty/")
intermediate_folder <- paste0(polarity, "_intermediate/")

falkor_metadata <- read.csv(file = "falkor_metadata.csv")
all_stans <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv")
found_stans <- read.csv(file = paste0(intermediate_folder, "found_stans.csv")) %>%
  left_join(all_stans, by=c(stan="Compound.Name")) %>%
  select(feature, stan, mzmed, rtmed, `old_stan`=Compound.Name_old)


addiso_peaks <- read.csv(paste0(pretty_folder, "addiso_peaks.csv")) %>%
  select(feature, mz, rt, into, file_name, M_area)
real_peaks <- read.csv(paste0(intermediate_folder, "real_peaks.csv")) %>%
  select(feature, mz, rt, into, file_name, M_area)
all_peaks <- rbind(addiso_peaks, real_peaks) %>%
  arrange("feature", "file_name")

quant <- read.csv("IS_HILIC-POS_Falkor.csv")
peaks_xcms <- found_stans %>%
  left_join(all_peaks, by="feature") %>%
  mutate(file_name=gsub(".mzML", "", .$file_name)) %>%
  left_join(quant, by=c(old_stan="Precursor.Ion.Name", file_name="Replicate.Name"))

ggplot(peaks_xcms) + 
  geom_point(aes(x=Area, y=into, color=stan, label=file_name)) +
  geom_abline(slope = 1) +
  facet_wrap(~old_stan, scales="free", ncol = 4) +
  theme(legend.position = "none") +
  xlab("Skyline area") +
  ylab("Custom XCMS area")
gp <- ggplot(peaks_xcms) + 
  geom_point(aes(x=log10(Area), y=log10(M_area), color=stan, label=file_name)) +
  geom_abline(slope = 1)
ggplotly(gp)




raw_data <- readRDS("stans_data.rds")
