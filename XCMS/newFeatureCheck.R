
library(tidyverse)
library(data.table)
library(gridExtra)

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

pmppm <- function(mass, ppm=4){c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))}



ms_files <- c("G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Blk_KM1906U14-Blk_C.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S62C1-25m_A.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S62C1-25m_B.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S62C1-25m_C.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S62C1-DCM_A.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S62C1-DCM_B.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S62C1-DCM_C.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S64C1-25m_A.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S64C1-25m_B.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S64C1-25m_C.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S64C1-DCM_A.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S64C1-DCM_B.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S64C1-DCM_C.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S77C1-25m_A.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S77C1-25m_B.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S77C1-25m_C.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S77C1-DCM_A.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S77C1-DCM_B.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S77C1-DCM_C.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S80C1-25m_A.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S80C1-25m_B.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S80C1-25m_C.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S80C1-DCM_A.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S80C1-DCM_B.mzML", 
              "G:\\My Drive\\FalkorFactor\\mzMLs\\190715_Smp_FK180310S80C1-DCM_C.mzML"
)

xdata_filled <- readRDS("XCMS/temp_data/current_xdata_filled.rds")

# raw_data <- lapply(ms_files, grabSingleFileData)
# raw_data <- lapply(seq_along(raw_data), function(x){
#   cbind(fileid=basename(ms_files)[x], raw_data[[x]])
# })
# raw_data <- do.call(rbind, raw_data)
# raw_data <- raw_data[raw_data$rt>60&raw_data$rt<1100,]
# saveRDS(raw_data, file = "MS1_data_frame")
MS1_data <- readRDS("XCMS/MS1_data_frame")
MS1_data$fileid <- as.character(MS1_data$fileid)
MS1_dt <- as.data.table(MS1_data)
rm(MS1_data)
split_MS1_dt <- split(MS1_dt, MS1_dt$fileid)

metadata <- data.frame(
  fileid=basename(ms_files),
  sample_group=c("Blank", "Sample")[c(1, rep(2, 24))],
  depth=c("Blank", "DCM", "25m")[c(1, rep(c(rep(2, 3), rep(3, 3)), 4))],
  spindir=c("Blank", "Cyclone", "Anticyclone")[c(1, rep(2, 12), rep(3, 12))],
  time=c("Blank", "Morning", "Afternoon")[c(1, rep(c(rep(2, 6), rep(3, 6)), 2))], stringsAsFactors = FALSE
)



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

interesting <- function(feature_data_i){
  missin_peaks <- feature_data_i %>%
    group_by(depth) %>%
    summarise(missin=sum(is.na(peak_area)))
  if(any(missin_peaks$missin>c(10, 10, 1))){
    return(NA)
  }
  # feature_data_i %>%
  #   aov(formula = peak_area~depth) %>%
  #   TukeyHSD() %>%
  #   `$`("depth") %>%
  #   `[`(,"p adj") %>%
  #   `<`(0.05) %>%
  #   any()
  t.test(feature_data_i$peak_area[feature_data_i$depth=="DCM"],
         feature_data_i$peak_area[feature_data_i$depth=="25m"]) %>%
    `$`("p.value") %>%
    `<`(0.01)
  # boxplot(cbind(feature_data_i$peak_area[feature_data_i$depth=="DCM"],
  #               feature_data_i$peak_area[feature_data_i$depth=="25m"]),
  #         names=c("DCM", "25m"), main=unique(feature_data_i$feature_num))
}

feature_data$interest <- feature_data %>%
  split(feature_data$feature_num) %>%
  sapply(interesting) %>% rep(each=25)
unique(feature_data$feature_num[feature_data$interest])

plotPeak <- function(feature_num_i){
  feature_info <- feature_data[feature_data$feature_num==feature_num_i,]
  feature_raw <- lapply(split(feature_info, feature_info$fileid), function(row_data){
    split_MS1_dt[[row_data$fileid]][
      rt%between%(c(row_data$rtmin, row_data$rtmax)+c(-50, 50)) &
        mz%between%c(row_data$mzmin, row_data$mzmax)]
  }) %>% do.call(what = rbind) %>%
    left_join(feature_info, by=c("fileid"))
  if(!nrow(feature_raw))return(ggplot()+ggtitle(label = feature_num_i)+theme_bw())
  
  gp <- ggplot(feature_raw) + geom_line(aes(x=rt, y=int, group=fileid, color=depth)) +
    scale_color_manual(values = c(`DCM`="#00FF0099", `25m`="#0000FF99", 
                                  `Blank`="#FF000099")) +
    geom_text(aes(x=Inf, y=Inf, label=unique(med_qscore)), hjust=1.5, vjust=1.5) +
    geom_text(aes(x=-Inf, y=Inf, label=round(mean(mz), digits = 5)), hjust=-0.5, vjust=1.5) +
    theme_bw() + theme(legend.position = "none") +
    ggtitle(label = feature_num_i)
  gp
}

pdf(file = "XCMS/newFeatureCheck.pdf", width = 17, height = 11)
for(i in seq(1, length(unique(feature_data$feature_num)), 16)){
  if(i+15>length(unique(feature_data$feature_num))){
    upperlim <- length(unique(feature_data$feature_num))
  } else {
    upperlim <- i+15
  }
  gps <- lapply(sprintf("FT%03d", i:upperlim), plotPeak)
  do.call("grid.arrange", gps)
}
dev.off()