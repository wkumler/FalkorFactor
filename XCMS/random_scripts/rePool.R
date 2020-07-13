library(dplyr)

col_names <- paste0(c("Mass [m/z]", "Formula [M]", "Formula type", "Species", "CS [z]", 
               "Polarity", "Start [min]", "End [min]", "Comment"), collapse = ",")

final_peaks <- read.csv("XCMS/pos_pretty/BMISed_peaks.csv")

features <- final_peaks %>%
  group_by(feature) %>%
  summarize(mzmed=median(mz), rtmin=min(rt)-60, rtmax=max(rt)+60) %>%
  as.data.frame()

csv_list <- lapply(split(features, features$feature), function(x){
  paste0(round(x$mzmed, digits = 4), ",,,,,Positive,", 
         floor(x$rtmin/60), ",", ceiling(x$rtmax/60), 
         ",", x$feature)
}) %>% unlist() %>% unname() %>% c(col_names, .)

writeLines(csv_list, con = "XCMS/random_scripts/fk_inclusion_list.csv")
