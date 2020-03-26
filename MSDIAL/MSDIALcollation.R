

#MSDIAL collation----
library(tidyverse)

MSDIAL_data <- lapply(c("pos", "neg"), function(pol){
  MSDIAL_outfiles <- list.files(paste0("MSDIAL/all_samples_", pol), 
                                pattern = ".txt", full.names = TRUE)
  MSDIAL_data <- lapply(MSDIAL_outfiles, function(filename){
    file_data <- suppressMessages(read_delim(filename, delim="\t", skip=4))
    data_long <- pivot_longer(file_data, cols = starts_with("190715_"),
                              values_to = filename %>% 
                                gsub(pattern = "^.*/", replacement = "") %>%
                                gsub(pattern = "_.*$", replacement = ""), 
                              names_to = "filename")
  })
  MSDIAL_pol <- Reduce(
    function(x, y, ...) suppressMessages(left_join(x, y, ...)), 
    MSDIAL_data
  )
  MSDIAL_pol$pol <- pol
  return(MSDIAL_pol)
}) %>% do.call(what = rbind) -> MSDIAL_data

MSDIAL_data <- MSDIAL_data %>% 
  group_by(`MS1 isotopic spectrum`) %>%
  summarise("S/N median" = median(SN)) %>% 
  left_join(MSDIAL_data, ., by='MS1 isotopic spectrum')
write.csv(x = MSDIAL_data, file = "MSDIAL/Collated_data.csv", row.names = FALSE)
