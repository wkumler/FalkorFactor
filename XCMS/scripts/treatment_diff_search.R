# Code that looks for differences in peak area across treatments
# Needs a peak list (typically produced by peakpicking.R)
# Calculates p-values between treatments
#  p-values are amended via the Holm-Bonferroni method to avoid multiple comparisons
#  Treatments are deduced from the file_name column of the peak list
# Data that only shows up in one treatment is given a p-value of 0.0000001

library(tidyverse)
final_peaks <- read.csv(file = "XCMS/data_pretty/final_peaks.csv")
final_features <- read.csv(file = "XCMS/data_pretty/final_features.csv")

diffreport <- function(peaks, pattern_a, pattern_b){
  split(peaks, peaks$feature) %>%
    lapply(function(x){
      DCM_areas <- x$M_area[grep(pattern = pattern_a, x$file_name)]
      m25_areas <- x$M_area[grep(pattern = pattern_b, x$file_name)]
      if(length(DCM_areas)<3|length(m25_areas)<3)return(c(0.0005, 1))
      c(pval=t.test(DCM_areas, m25_areas)$p.value, 
        diff=mean(DCM_areas)/mean(m25_areas))
    }) %>% do.call(what = rbind) %>% as.data.frame() %>%
    mutate(feature=rownames(.)) %>%
    arrange(pval) %>%
    mutate(p_adj=p.adjust(pval)) %>%
    filter(p_adj<0.05) %>%
    arrange(p_adj) %>%
    mutate(which_enriched=c(pattern_a, pattern_b)[as.numeric(diff>1)+1]) %>%
    split(.$which_enriched) %>%
    lapply(function(x){x[names(x)!="which_enriched"]})
}

# Depth data!
diffreport(final_peaks, pattern_a = "DCM", pattern_b = "25m")

# Diel data!
diffreport(final_peaks, pattern_a = "62|77", pattern_b = "64|80")

# Direction data!
diffreport(final_peaks, pattern_a = "62|64", pattern_b = "77|80")



left_join(diffreport(final_peaks, pattern_a = "DCM", pattern_b = "25m")[[1]],
      final_features, by="feature") %>% arrange(p_adj)

filter(final_peaks, feature=="FT305")


feature_formulas <- read.csv(file = "XCMS/data_pretty/feature_formulas.csv")
