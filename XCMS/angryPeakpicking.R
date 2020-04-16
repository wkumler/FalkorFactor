# Angry peakpicking
library(xcms)
library(tidyverse)
library(data.table)

register(BPPARAM = SerialParam())

ms_files <- "mzMLs" %>%
  list.files(pattern = ".mzML", full.names = TRUE) %>%
  normalizePath() %>%
  `[`(!grepl("Fullneg|Fullpos|QC-KM1906", x = .))
metadata <- data.frame(
  fileid=basename(ms_files),
  sample_group=c("Blank", "Pooled", "Sample", "Std")[c(1, rep(2, 6), rep(3, 24), rep(4, 10))],
  depth=c("Blank", "Pooled", "DCM", "25m", "Std")[c(1, rep(2, 6), rep(c(rep(3, 3), rep(4, 3)), 4), rep(5, 10))],
  spindir=c("Blank", "Pooled", "Cyclone", "Anticyclone", "Std")[c(1, rep(2, 6), rep(3, 12), rep(4, 12), rep(5, 10))],
  time=c("Blank", "Pooled", "Morning", "Afternoon", "Std")[c(1, rep(2, 6), rep(c(rep(3, 6), rep(4, 6)), 2), rep(5, 10))]
)
raw_data <- readMSData(files = ms_files, pdata = metadata, 
                       mode = "onDisk", verbose = TRUE)

all_data <- as.data.table(readRDS(file = "XCMS/MS1_data_frame"))

given_mass <- 143.082053
given_mass <- 118.0865
given_rt <- 450
given_rt <- 350

rtr <- (given_rt+c(-100, 100))
mzr <- pmppm(given_mass, 10)
group_colors <- c("red", "black", "steelblue", NA)
chr_raw <- xcms::chromatogram(raw_data, mz = mzr, rt = rtr)
plot(chr_raw, col = group_colors[chr_raw$sample_group])

cwp <- CentWaveParam(ppm = 5, peakwidth = c(10, 20), 
                     snthresh = 0, prefilter = c(0,0), 
                     integrate = 2, mzCenterFun = "wMean", 
                     mzdiff = 0.0001, fitgauss = FALSE, 
                     noise = 0, firstBaselineCheck = FALSE, 
                     extendLengthMSW = TRUE)
xchr <- findChromPeaks(chr_raw, param = cwp)
chrom_peaks <- as.data.frame(chromPeaks(xchr)) %>% 
  mutate(fileid=column) %>% filter(fileid<26)

eics <- all_data[mz%between%pmppm(given_mass, 5)&
                   rt%between%(given_rt+c(-100, 100))]
ggplot() + 
  geom_line(data = eics, aes(x=rt, y=int, group=fileid)) + 
  geom_vline(data = chrom_peaks, aes(xintercept=rtmin+1), col="red") +
  geom_vline(data = chrom_peaks, aes(xintercept=rtmax-1), col="red") #+
  facet_wrap(~fileid, scales = "free_y")
chrom_peaks <- mutate(chrom_peaks, mzmin=min(eics$mz), mzmax=max(eics$mz))

lapply(seq_len(nrow(chrom_peaks)), FUN=xcmsQscoreCalculator, 
       xcms_peakdf=chrom_peaks, file_data_table=eics, 
       qscoreCalculator=qscoreCalculator) %>%
  do.call(what = rbind) %>%
  cbind(chrom_peaks, .) %>%
  filter(fileid==5) %>% arrange(rt)

head(xcms::rtime(xchr[[5]]))
head(eics[fileid==5])
