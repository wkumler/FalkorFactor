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
) %>% new(Class = "NAnnotatedDataFrame")

# Functions ----
pmppm <- function(mass, ppm=4){c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))}
qscoreCalculator <- function(eic){
  #Check for bogus EICs
  if(nrow(eic)<5){
    return(data.frame(SNR=0, peak_cor=0, qscore=0, skew=0))
  }
  #Calculate where each rt would fall on a beta distribution (accounts for missed scans)
  scaled_rts <- (eic$rt-min(eic$rt))/(max(eic$rt)-min(eic$rt))
  
  # Create a couple different skews and test fit
  possible_skews <- c(2,2.5,3,4,5,7)
  best_skew <- possible_skews[which.max(sapply(possible_skews, function(x){
    cor(dbeta(scaled_rts, shape1 = x, shape2 = 5), eic$int)
  }))]
  perf_peak <- dbeta(scaled_rts, shape1 = best_skew, shape2 = 5)
  peak_cor <- cor(perf_peak, eic$int)
  
  #Calculate the normalized residuals
  residuals <- eic$int/max(eic$int)-perf_peak/max(perf_peak)
  #Calculate the minimum SD, after normalizing for any shape discrepancy
  old_res_sd <- sd(residuals)
  norm_residuals <- diff(residuals)
  new_res_sd <- sd(norm_residuals)
  while(new_res_sd<old_res_sd){
    old_res_sd <- new_res_sd
    norm_residuals <- diff(residuals)
    new_res_sd <- sd(residuals)
  }
  #Calculate SNR
  SNR <- (max(eic$int)-min(eic$int))/sd(norm_residuals*max(eic$int))
  #Return the quality score
  output <- data.frame(SNR, peak_cor, 
                       qscore=SNR*peak_cor^4*log10(max(eic$int)), 
                       skew=best_skew)
  return(output)
}
xcmsQscoreCalculator <- function(df_row, xcms_peakdf, all_data_table, 
                                 qscoreCalculator = qscoreCalculator){
  #Extract the relevant EIC
  peak_row_data <- xcms_peakdf[df_row, ]
  file_data_table <- all_data_table[all_data_table$fileid==xcms_peakdf[
    df_row, "fileid"]]
  eic <- file_data_table[rt %between% c(peak_row_data$rtmin, peak_row_data$rtmax)&
                           mz %between% c(peak_row_data$mzmin, peak_row_data$mzmax)]
  return(qscoreCalculator(eic))
}


# Angriness ----
raw_data <- readMSData(files = ms_files, pdata = metadata, 
                       mode = "onDisk", verbose = TRUE)

all_data <- as.data.table(readRDS(file = "XCMS/MS1_data_frame"))



given_mass <- 118.0865
given_rt <- 350
rtrbet <- (given_rt+c(-100, 100))
mzrbet <- pmppm(given_mass, 5)
eicsbet <- all_data[mz%between%mzrbet & rt%between%rtrbet]
group_colors <- c("red", "black", "steelblue", NA)
chr_bet <- xcms::chromatogram(raw_data, mz = mzrbet, rt = rtrbet)
plot(chr_bet, col = group_colors[chr_bet$sample_group])

given_mass <- 143.082053
given_rt <- 450
rtrect <- (given_rt+c(-100, 100))
mzrect <- pmppm(given_mass, 5)
eicsect <- all_data[mz%between%mzrect & rt%between%rtrect]
group_colors <- c("red", "black", "steelblue", NA)
chr_ect <- xcms::chromatogram(raw_data, mz = mzrect, rt = rtrect)
plot(chr_ect, col = group_colors[chr_ect$sample_group])

given_mass <- 123.040504455566
given_rt <- 590
rtrbleh <- (given_rt+c(-100, 100))
mzrbleh <- pmppm(given_mass, 5)
eicsbleh <- all_data[mz%between%mzrbleh & rt%between%rtrbleh]
group_colors <- c("red", "black", "steelblue", NA)
chr_bleh <- xcms::chromatogram(raw_data, mz = mzrbleh, rt = rtrbleh)
plot(chr_bleh, col = group_colors[chr_bleh$sample_group])



cwp <- CentWaveParam(ppm = 5, peakwidth = c(15, 15), 
                     snthresh = 0, prefilter = c(0,0), 
                     integrate = 2, mzCenterFun = "wMean", 
                     mzdiff = 0.001, fitgauss = FALSE, 
                     noise = 0, firstBaselineCheck = FALSE, 
                     extendLengthMSW = TRUE)

xchr <- findChromPeaks(chr_bet, param = cwp)
chrom_peaks <- as.data.frame(chromPeaks(xchr), stringsAsFactors = FALSE) %>%
  mutate(fileid=metadata@data$fileid[column]) %>%
  mutate(fileid=as.character(fileid)) %>%
  filter(fileid%in%unique(eicsbet$fileid)) %>%
  mutate(mzmin=min(eicsbet$mz), mzmax=max(eicsbet$mz))

xcmsQscoreCalculator(df_row = 3, xcms_peakdf = chrom_peaks, 
                     all_data_table = eicsbet, 
                     qscoreCalculator = qscoreCalculator)


qscore_peaks <- lapply(seq_len(nrow(chrom_peaks)), FUN=xcmsQscoreCalculator,
       xcms_peakdf=chrom_peaks, all_data_table=eicsbet,
       qscoreCalculator=qscoreCalculator) %>%
  do.call(what = rbind) %>%
  cbind(chrom_peaks, .) %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  select("rt", "sn", "fileid", "SNR", "peak_cor", "qscore")

ggplot() + 
  geom_line(data = eicsbet, aes(x=rt, y=int, group=fileid)) + 
  facet_wrap(~fileid, scales = "free_y") +
  geom_vline(data = chrom_peaks, aes(xintercept=rtmin+1), col="red") +
  geom_vline(data = chrom_peaks, aes(xintercept=rtmax-1), col="red") +
  geom_text(data = qscore_peaks, aes(x=rt, y=0, label=round(qscore)))
ggsave(filename = "~/../Desktop/RplotBet.pdf", device = "pdf", 
       width = 17, height=11)


xchr <- findChromPeaks(chr_ect, param = cwp)
chrom_peaks <- as.data.frame(chromPeaks(xchr), stringsAsFactors = FALSE) %>%
  mutate(fileid=metadata@data$fileid[column]) %>%
  mutate(fileid=as.character(fileid)) %>%
  filter(fileid%in%unique(eicsect$fileid)) %>%
  mutate(mzmin=min(eicsect$mz), mzmax=max(eicsect$mz))

xcmsQscoreCalculator(df_row = 3, xcms_peakdf = chrom_peaks, 
                     all_data_table = eicsect, 
                     qscoreCalculator = qscoreCalculator)


qscore_peaks <- lapply(seq_len(nrow(chrom_peaks)), FUN=xcmsQscoreCalculator,
                       xcms_peakdf=chrom_peaks, all_data_table=eicsect,
                       qscoreCalculator=qscoreCalculator) %>%
  do.call(what = rbind) %>%
  cbind(chrom_peaks, .) %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  select("rt", "sn", "fileid", "SNR", "peak_cor", "qscore")

ggplot() + 
  geom_line(data = eicsect, aes(x=rt, y=int, group=fileid)) + 
  facet_wrap(~fileid, scales = "free_y") +
  geom_vline(data = chrom_peaks, aes(xintercept=rtmin+1), col="red") +
  geom_vline(data = chrom_peaks, aes(xintercept=rtmax-1), col="red") +
  geom_text(data = qscore_peaks, aes(x=rt, y=0, label=round(qscore)))
ggsave(filename = "~/../Desktop/RplotEct.pdf", device = "pdf", 
       width = 17, height=11)


xchr <- findChromPeaks(chr_bleh, param = cwp)
chrom_peaks <- as.data.frame(chromPeaks(xchr), stringsAsFactors = FALSE) %>%
  mutate(fileid=metadata@data$fileid[column]) %>%
  mutate(fileid=as.character(fileid)) %>%
  filter(fileid%in%unique(eicsbleh$fileid)) %>%
  mutate(mzmin=min(eicsbleh$mz), mzmax=max(eicsbleh$mz))

xcmsQscoreCalculator(df_row = 3, xcms_peakdf = chrom_peaks, 
                     all_data_table = eicsbleh, 
                     qscoreCalculator = qscoreCalculator)


qscore_peaks <- lapply(seq_len(nrow(chrom_peaks)), FUN=xcmsQscoreCalculator,
                       xcms_peakdf=chrom_peaks, all_data_table=eicsbleh,
                       qscoreCalculator=qscoreCalculator) %>%
  do.call(what = rbind) %>%
  cbind(chrom_peaks, .) %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  select("rt", "sn", "fileid", "SNR", "peak_cor", "qscore")

ggplot() + 
  geom_line(data = eicsbleh, aes(x=rt, y=int, group=fileid)) + 
  facet_wrap(~fileid, scales = "free_y") +
  geom_vline(data = chrom_peaks, aes(xintercept=rtmin+1), col="red") +
  geom_vline(data = chrom_peaks, aes(xintercept=rtmax-1), col="red") +
  geom_text(data = qscore_peaks, aes(x=rt, y=0, label=round(qscore)))
ggsave(filename = "~/../Desktop/RplotBleh.pdf", device = "pdf", 
       width = 17, height=11)
