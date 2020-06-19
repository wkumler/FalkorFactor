# Code that looks for differences in peak area across treatments
# Needs a peak list (typically produced by peakpicking.R)
# Calculates p-values between treatments
#  p-values are amended via the Holm-Bonferroni method to avoid multiple comparisons
#  Treatments are deduced from the file_name column of the peak list
# Data that only shows up in one treatment is given a p-value of 0

library(tidyverse)
final_peaks <- read.csv(file = "XCMS/data_pretty/final_peaks.csv")
final_features <- read.csv(file = "XCMS/data_pretty/final_features.csv")
feature_formulas <- read.csv(file = "XCMS/data_pretty/feature_formulas.csv")
feature_putatives <- "XCMS/data_intermediate/sirius_project/output_dir/" %>%
  paste0("compound_identifications.tsv") %>%
  readLines() %>%
  strsplit(split = "\t")
header <- feature_putatives[[1]]
feature_putatives <- feature_putatives[-1] %>%
  do.call(what = rbind) %>%
  as.data.frame() %>%
  `names<-`(header) %>%
  select(id, Confidence_Score, name) %>%
  mutate(feature=gsub(pattern = "^[[:digit:]]+_", "", .$id)) %>%
  mutate(feature=gsub(pattern = "_FEATURE1$", "", .$feature)) %>%
  mutate(Confidence_Score=round(as.numeric(Confidence_Score), digits = 3)) %>%
  select(-id)

diffreport <- function(peaks, pattern_a, pattern_b){
  split(peaks, peaks$feature) %>%
    lapply(function(x){
      DCM_areas <- x$BMISed_area[grep(pattern = pattern_a, x$file_name)]
      m25_areas <- x$BMISed_area[grep(pattern = pattern_b, x$file_name)]
      if(length(DCM_areas)<3|length(m25_areas)<3)return(c(0, 1))
      c(pval=t.test(DCM_areas, m25_areas)$p.value, 
        diff=mean(DCM_areas)/mean(m25_areas))
    }) %>% do.call(what = rbind) %>% as.data.frame() %>%
    mutate(feature=rownames(.)) %>%
    arrange(pval) %>%
    mutate(p_adj=p.adjust(pval, method = "BH")) %>%
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



mega_df <- final_peaks %>%
  diffreport(pattern_a = "DCM", pattern_b = "25m") %>%
  #diffreport(pattern_a = "62|77", pattern_b = "64|80") %>%
  #diffreport(pattern_a = "62|64", pattern_b = "77|80") %>%
  `[[`(1) %>%
  left_join(final_features, by="feature") %>% 
  arrange(p_adj) %>%
  left_join(feature_formulas, by=c("feature", "mzmed", "rtmed", "avgarea")) %>%
  select(feature, p_adj, mzmed, rtmed, avgarea, C13, S34, N15, O18, formula) %>%
  left_join(feature_putatives, by="feature")
mega_df <- final_peaks %>%
  diffreport(pattern_a = "DCM", pattern_b = "25m") %>%
  do.call(what = rbind) %>%
  left_join(final_features, by="feature") %>%
  arrange(p_adj) %>%
  select(feature, p_adj, mzmed, rtmed, BMIS_avg)
as_tibble(mega_df)

