
# Setup things ----
library(tidyverse)
library(data.table)
library(pbapply)


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



# Load raw data ----
ms_files <- "mzMLs" %>%
  list.files(pattern = ".mzML", full.names = TRUE) %>%
  normalizePath() %>%
  `[`(!grepl("Fullneg|Fullpos|QC-KM1906", x = .))

dt <- pblapply(seq_along(ms_files), function(x){
  v <- grabSingleFileData(ms_files[x]) %>%
    cbind(file_name=basename(ms_files[x]), .)
  return(v)
}) %>% do.call(what = rbind) %>% as.data.table()

metadata <- data.frame(
  fileid=1:41, 
  filenames=basename(ms_files),
  sample_group=c("Blank", "Pooled", "Sample", "Std")[c(1, rep(2, 6), rep(3, 24), rep(4, 10))],
  depth=c("Blank", "Pooled", "DCM", "25m", "Std")[c(1, rep(2, 6), rep(c(rep(3, 3), rep(4, 3)), 4), rep(5, 10))],
  spindir=c("Blank", "Pooled", "Cyclone", "Anticyclone", "Std")[c(1, rep(2, 6), rep(3, 12), rep(4, 12), rep(5, 10))],
  time=c("Blank", "Pooled", "Morning", "Afternoon", "Std")[c(1, rep(2, 6), rep(c(rep(3, 6), rep(4, 6)), 2), rep(5, 10))]
)



# Load MSDIAL data ----
features <- read_delim(file = "MSDIAL/all_samples_pos/Area_0_20202121116.txt", 
                       delim = "\t", skip = 4)
clean_features <- features %>% 
  arrange(desc(`S/N average`)) %>% 
  
  select(starts_with("Average"), "S/N average")

i <- 1
while(TRUE){
  mzr <- pmppm(clean_features[i,"Average Mz"])
  rtr <- clean_features[i,"Average Rt(min)"] %>% c(.-1, .+1) %>% `[`(-1) %>% 
    unlist() %>% `*`(60)
  eic <- dt[rt%between%rtr&mz%between%mzr]
  anno_eic <- left_join(eic, metadata, by=c(file_name="filenames"))
  gp <- ggplot(anno_eic) + geom_line(aes(x=rt, y=int, group=file_name, color=sample_group)) +
    ggtitle(paste("m/z:", clean_features[i,"Average Mz"], 
                  "   s/n:", clean_features[i,"S/N average"])) +
    scale_color_manual(values = c(Blank="black", Pooled="red", Sample="blue", Std="grey"))
  print(gp)
  readline(prompt = "Press to continue")
  i <- i+1
}
