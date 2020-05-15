# Code that looks for differences in peak area across treatments
# Needs a peak list (typically produced by peakpicking.R)
# Calculates p-values between treatments
#  Treatments are deduced from the file_name column of the peak list
# Data that only shows up in one treatment is given a p-value of 0.0005



library(tidyverse)
final_peaks <- read.csv(file = "XCMS/data_pretty/final_peaks.csv")



# Depth data!
final_diffreport <- split(final_peaks, final_peaks$feature) %>%
  lapply(function(x){
    DCM_areas <- x$M_area[grep(pattern = "DCM", x$file_name)]
    m25_areas <- x$M_area[grep(pattern = "25m", x$file_name)]
    if(length(DCM_areas)<3|length(m25_areas)<3)return(c(0.0005, 1))
    c(pval=t.test(DCM_areas, m25_areas)$p.value, 
      diff=mean(DCM_areas)/mean(m25_areas))
  }) %>% do.call(what = rbind) %>% as.data.frame() %>%
  mutate(feature=rownames(.)) %>% 
  group_by(feature) %>%
  summarize(mzmed=median(mz), rtmed=median(rt), avgarea=mean(area)) %>%
  arrange(mzmed)
DCM_enriched <- final_diffreport %>%
  filter(diff>1) %>%
  filter(pval<0.01) %>%
  arrange(pval)
surface_enriched <- final_diffreport %>%
  filter(diff<1) %>%
  filter(pval<0.01) %>%
  arrange(pval)



# Diel data!
final_diffreport <- split(final_peaks, final_peaks$feature) %>%
  lapply(function(x){
    AM_areas <- x$M_area[grep(pattern = "62|77", x$file_name)]
    PM_areas <- x$M_area[grep(pattern = "64|80", x$file_name)]
    if(length(AM_areas)<3|length(PM_areas)<3)return(c(0.0005, 1))
    c(pval=t.test(AM_areas, PM_areas)$p.value, 
      diff=mean(AM_areas)/mean(PM_areas))
  }) %>% do.call(what = rbind) %>% as.data.frame() %>%
  mutate(feature=rownames(.)) %>% 
  group_by(feature) %>%
  summarize(mzmed=median(mz), rtmed=median(rt), avgarea=mean(area)) %>%
  arrange(mzmed)
AM_enriched <- final_diffreport %>%
  filter(diff>1) %>%
  filter(pval<0.01) %>%
  arrange(pval)
PM_enriched <- final_diffreport %>%
  filter(diff<1) %>%
  filter(pval<0.01) %>%
  arrange(pval)
