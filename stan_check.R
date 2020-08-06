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
pretty_folder <- paste0("XCMS/", polarity, "_pretty/")
intermediate_folder <- paste0("XCMS/", polarity, "_intermediate/")

falkor_metadata <- read.csv(file = "XCMS/falkor_metadata.csv")
found_stans <- read.csv(file = paste0(intermediate_folder, "found_stans.csv"))
all_stans <- read.csv("https://raw.githubusercontent.com/IngallsLabUW/Ingalls_Standards/master/Ingalls_Lab_Standards_NEW.csv")
found_stans <- found_stans %>%
  left_join(all_stans, by=c(stan="Compound.Name")) %>%
  select(names(found_stans), old_stan=Compound.Name_old)


addiso_peaks <- read.csv(paste0(pretty_folder, "addiso_peaks.csv")) %>%
  select(feature, mz, rt, into, file_name, M_area)
real_peaks <- read.csv(paste0(intermediate_folder, "real_peaks.csv")) %>%
  select(feature, mz, rt, into, file_name, M_area)
all_peaks <- rbind(addiso_peaks, real_peaks) %>%
  arrange("feature", "file_name")
found_stans <- found_stans %>%
  left_join(all_peaks) %>%
  select(old_stan, file_name, into, M_area)

quant <- read.csv("C:/Users/willi/Downloads/IS_HILIC-POS_Falkor.csv")
found_stans %>%
  mutate(file_name=gsub(".mzML", "", .$file_name)) %>%
  left_join(quant, by=c(old_stan="Precursor.Ion.Name", file_name="Replicate.Name"))







arseno_xcms <- filter(all_peaks, feature=="FT419") %>%
  select(file_name, into, M_area) %>%
  mutate(file_name=gsub(".mzML", "", .$file_name)) %>%
  arrange(file_name)
arseno_quant <- filter(quant, Precursor.Ion.Name=="Arsenobetaine, 13C2") %>%
  select(file_name=Replicate.Name, Area) %>%
  arrange(file_name)

v <- left_join(arseno_quant, arseno_xcms) #%>%
  ggplot() +
  geom_point(aes(x=Area, y=M_area)) +
  geom_abline(slope = 1) +
  coord_fixed()
