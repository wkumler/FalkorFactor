# Code to analyze the output of MSDIAL's peakpicking when applied to Falkor's
# pooled samples to get a sense of which peaks are real

# File source: PooledOnly_Height.txt
# From: MSDIAL PooledOnly.mtd2 -> export -> alignment -> Height, txt, centroid

# Goal: Output a file of peak m/z and rt to check out peak shape in XCMS


library(xcms)
library(tidyverse)

will_plot <- function(MSnExp_obj){
  x <- as(MSnExp_obj, "data.frame")
  x <- split(x, x$file)
  lapply(x, function(y){
    bpi <- unlist(lapply(split(y$i, y$rt), max, na.rm = TRUE))
    hclcols = hcl.colors(20)[cut(bpi, breaks = 20)]
    plot(as.numeric(names(bpi)), bpi, col=hclcols, pch=19)
  })
}

all_peaks <- read.table("MSDIAL/Processing/Height_0_2019104109.txt", skip = 4, 
                        sep = "\t", header = T)[,c(1,2,3,4,29:56)]



diff_vals <- numeric(nrow(all_peaks))
dcm_cols <- grep("DCM",names( all_peaks))
m25_cols <- grep("25m",names( all_peaks))
for(i in seq_len(nrow(all_peaks))){
  diff_vals[i] <- t.test(all_peaks[i, dcm_cols], all_peaks[i, m25_cols])$p.value
}
all_peaks$DCM_25m_pvals <- diff_vals


diff_vals <- numeric(nrow(all_peaks))
ccw_cols <- grep("62|64", names(all_peaks))
cw_cols <- grep("77|80", names(all_peaks))
for(i in seq_len(nrow(all_peaks))){
  diff_vals[i] <- t.test(all_peaks[i, cw_cols], all_peaks[i, ccw_cols])$p.value
}
all_peaks$cw_ccw_pvals <- diff_vals

save(all_peaks, file = "Preprocessing_all_peaks")
load("Preprocessing_all_peaks")

# for(row in which(all_peaks$DCM_25m_pvals<1e-04)){
#   boxplot(as.numeric(all_peaks[row, dcm_cols]), as.numeric(all_peaks[row, m25_cols]))
# }
# for(row in which(all_peaks$cw_ccw_pvals<1e-03)){
#   boxplot(as.numeric(all_peaks[row, ccw_cols]), as.numeric(all_peaks[row, cw_cols]))
# }

cw_ccw_enriched <- subset(all_peaks, cw_ccw_pvals<0.001)




# Check interesting peak shapes

mzml_files <- list.files("mzMLs", full.names = T)[c(17:50)]
raw_data <- readMSData(files = mzml_files, msLevel. = 1, mode = "onDisk")
par(mfrow=c(2,2))
par(mar=c(2.1, 2.1, 0.1, 0.1))
for(i in 1:10){
  rt_i <- known_peaks$Average.Rt.min.[i]*60
  mz_i <- known_peaks$Average.Mz[i]
  raw_data %>%
    filterMz(mz_i+c(mz_i*-2.5, mz_i*2.5)/1000000) %>%
    filterRt(c(rt_i-100, rt_i+100)) %>%
    will_plot()
}
