# Angry peakpicking

xdata_cor <- readRDS(file = "XCMS/temp_data/current_xdata_cor.rds")
all_data <- as.data.table(readRDS(file = "MS1_data_frame"))
peakdf_qscored <- as.matrix(read.csv(file = "XCMS/temp_data/peakdf_qscored.csv"))



xdata_cor %>% featureDefinitions() %>% as.data.frame() %>% 
  filter(mzmed%between%pmppm(90.05550455, 5))

chromPeaks(xdata_cor) %>% as.data.frame() %>% 
  filter(mz%between%pmppm(90.05550455, 10)) %>%
  arrange(rt)

peakdf_qscored %>% as.data.frame() %>% 
  filter(mz%between%pmppm(90.05550455, 10)) %>%
  ggplot() + geom_rect(aes(xmin=rtmin, xmax=rtmax, ymin=mzmin, 
                           ymax=mzmax, fill=log10(maxo)), color="black") +
  facet_wrap(~sample) + scale_fill_viridis_c() +
  xlim(c(500, 650))

eics <- all_data[mz%between%pmppm(90.05550455, 10)&
                   rt%between%c(500,650)]
ggplot(eics) + geom_line(aes(x=rt, y=int)) + 
  facet_wrap(~fileid, scales = "free_y")
