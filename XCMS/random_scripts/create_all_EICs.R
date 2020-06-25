

# Setup things ----
library(tidyverse)
library(data.table)
library(gridExtra)
library(pbapply)

#Only for RScript running
setwd(r"(G:\My Drive\FalkorFactor)")

ms_files <- "mzMLs" %>%
  list.files(pattern = "Blk|Smp", full.names = TRUE) %>%
  normalizePath()
xdata_filled <- readRDS("XCMS/data_intermediate/current_xdata_filled.rds")



# Functions ----
pmppm <- function(mass, ppm=4){c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))}
grabSingleFileData <- function(filename){
  msdata <- mzR:::openMSfile(filename)
  fullhd <- mzR::header(msdata)
  spectra_list <- lapply(seq_len(nrow(fullhd)), function(x){
    given_peaks <- mzR::peaks(msdata, x)
    rtime <- fullhd[x, "retentionTime"]
    return(cbind(rtime, given_peaks))
  })
  all_data <- `names<-`(as.data.frame(do.call(rbind, spectra_list)), 
                        c("rt", "mz", "int"))
  return(all_data)
}
plotPeak <- function(feature_num_i){
  feature_info <- feature_data[feature_data$feature_num==feature_num_i,]
  feature_raw <- lapply(split(feature_info, feature_info$fileid), function(row_data){
    split_MS1_dt[[row_data$fileid]][
      rt%between%(c(row_data$rtmin, row_data$rtmax)+c(-30, 80)) &
        mz%between%c(row_data$mzmin, row_data$mzmax)]
  }) %>% do.call(what = rbind) %>%
    left_join(feature_info, by=c("fileid"))
  if(!nrow(feature_raw))return(ggplot()+ggtitle(label = feature_num_i)+theme_bw())
  
  gp <- ggplot(feature_raw) + geom_line(aes(x=rt, y=int, group=fileid, color=depth)) +
    scale_color_manual(values = c(`DCM`="#00FF0099", `25m`="#0000FF99", 
                                  `Blank`="#FF000099", `Standard`="#00000099", 
                                  `Pooled`="#00FFFF99")) +
    geom_text(aes(x=Inf, y=Inf, label=unique(med_qscore)), hjust=1.5, vjust=1.5) +
    geom_text(aes(x=-Inf, y=Inf, label=round(mean(mz), digits = 5)), hjust=-0.5, vjust=1.5) +
    theme_bw() + theme(legend.position = "none") +
    ggtitle(label = feature_num_i)
  gp
}



# Read in raw data ----
raw_data <- pblapply(ms_files, grabSingleFileData)
raw_data <- pblapply(seq_along(raw_data), function(x){
  cbind(fileid=basename(ms_files)[x], raw_data[[x]])
}) %>% do.call(what = rbind) %>%
  `[`(.$rt>60&.$rt<1100,)
saveRDS(raw_data, file = "XCMS/data_intermediate/MS1_data_frame.rds")
MS1_dt <- readRDS("XCMS/data_intermediate/MS1_data_frame.rds") %>%
  mutate(fileid=as.character(.$fileid)) %>%
  as.data.table()
split_MS1_dt <- split(MS1_dt, MS1_dt$fileid)

metadata <- data.frame(
  fileid=basename(ms_files),
  sample_group=c("Blank", "Sample")[c(1, rep(2, 24))],
  depth=c("Blank", "DCM", "25m")[c(1, rep(c(rep(2, 3), rep(3, 3)), 4))],
  spindir=c("Blank", "Cyclone", "Anticyclone")[c(1, rep(2, 12), rep(3, 12))],
  time=c("Blank", "Morning", "Afternoon")[c(1, rep(c(rep(2, 6), rep(3, 6)), 2))], stringsAsFactors = FALSE
)



# Calculate quality scores ----
feature_qscores <- xdata_filled %>%
  xcms::featureDefinitions() %>%
  `$`("peakidx") %>%
  sapply(FUN = function(x){
    median(xcms::chromPeaks(xdata_filled)[x,"sn"], na.rm = TRUE)
  }) %>%
  data.frame(feature_num=rownames(xcms::featureDefinitions(xdata_filled)),
             med_qscore=., stringsAsFactors = FALSE)

feature_data <- xdata_filled %>% 
  xcms::featureValues() %>% 
  cbind(feature_num=rownames(.), .) %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  pivot_longer(cols = -feature_num, 
               names_to = "fileid", values_to = "peak_area") %>%
  mutate(peak_area=as.numeric(peak_area)) %>%
  left_join(metadata %>% mutate(fileid=as.character(fileid)), by="fileid") %>%
  filter(sample_group%in%c("Blank", "Sample")) %>%
  left_join(xdata_filled %>% 
              xcms::featureDefinitions() %>% 
              as.data.frame(stringsAsFactors=FALSE) %>%
              select(mzmin,mzmax,rtmin,rtmax) %>% 
              cbind(feature_num=rownames(.), .) %>%
              mutate(feature_num=as.character(feature_num)), 
            by="feature_num") %>%
  left_join(feature_qscores, by="feature_num")



# Plot each chromatogram and save to pdf ----
pdf(file = "XCMS/data_pretty/all_feature_EICs.pdf", width = 17, height = 11)
pb <- txtProgressBar(min = 0, max = length(unique(feature_data$feature_num)), style = 3)
for(i in seq(1, length(unique(feature_data$feature_num)), 16)){
  if(i+15>length(unique(feature_data$feature_num))){
    upperlim <- length(unique(feature_data$feature_num))
  } else {
    upperlim <- i+15
  }
  gps <- lapply(sprintf("FT%03d", i:upperlim), plotPeak)
  do.call("grid.arrange", gps)
  setTxtProgressBar(pb, i)
}
close(pb)
dev.off()



# Repeat for all files and RT's ----
ms_files <- "mzMLs" %>%
  list.files(pattern = ".mzML", full.names = TRUE) %>%
  setdiff(list.files(path = "mzMLs", pattern = "neg|pos", full.names = TRUE)) %>%
  normalizePath()
raw_data <- pblapply(ms_files, grabSingleFileData)
raw_data <- pblapply(seq_along(raw_data), function(x){
  cbind(fileid=basename(ms_files)[x], raw_data[[x]])
}) %>% do.call(what = rbind)
saveRDS(raw_data, file = "XCMS/data_intermediate/MS1_data_frame_all.rds")
MS1_dt <- readRDS("XCMS/data_intermediate/MS1_data_frame_all.rds") %>%
  mutate(fileid=as.character(.$fileid)) %>%
  as.data.table()
split_MS1_dt <- split(MS1_dt, MS1_dt$fileid)
metadata <- data.frame(
  fileid=basename(ms_files),
  sample_group=c("Blank", "QC", "Pooled", "Sample", "Standard")[
    c(1, rep(2, 3), rep(3, 6), rep(4, 24), rep(5, 10))],
  depth=c("Blank", "QC", "Pooled", "25m", "DCM", "Standard")[
    c(1, rep(2, 3), rep(3, 6), rep(c(rep(4, 3), rep(5, 3)), 4), rep(6, 10))])
feature_qscores <- xdata_filled %>%
  xcms::featureDefinitions() %>%
  `$`("peakidx") %>%
  sapply(FUN = function(x){
    median(xcms::chromPeaks(xdata_filled)[x,"sn"], na.rm = TRUE)
  }) %>%
  data.frame(feature_num=rownames(xcms::featureDefinitions(xdata_filled)),
             med_qscore=., stringsAsFactors = FALSE)
feature_data <- xdata_filled %>% 
  xcms::featureValues() %>% 
  cbind(feature_num=rownames(.), .) %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  pivot_longer(cols = -feature_num, 
               names_to = "fileid", values_to = "peak_area") %>%
  mutate(peak_area=as.numeric(peak_area)) %>%
  left_join(metadata %>% mutate(fileid=as.character(fileid)), by="fileid") %>%
  left_join(xdata_filled %>% 
              xcms::featureDefinitions() %>% 
              as.data.frame(stringsAsFactors=FALSE) %>%
              select(mzmin,mzmax,rtmin,rtmax) %>% 
              cbind(feature_num=rownames(.), .) %>%
              mutate(feature_num=as.character(feature_num)), 
            by="feature_num") %>%
  left_join(feature_qscores, by="feature_num")
pdf(file = "XCMS/data_pretty/all_feature_EICs_w_stans.pdf", width = 17, height = 11)
pb <- txtProgressBar(min = 0, max = length(unique(feature_data$feature_num)), style = 3)
for(i in seq(1, length(unique(feature_data$feature_num)), 16)){
  if(i+15>length(unique(feature_data$feature_num))){
    upperlim <- length(unique(feature_data$feature_num))
  } else {
    upperlim <- i+15
  }
  gps <- lapply(sprintf("FT%03d", i:upperlim), plotPeak)
  do.call("grid.arrange", gps)
  setTxtProgressBar(pb, i)
}
close(pb)
dev.off()
