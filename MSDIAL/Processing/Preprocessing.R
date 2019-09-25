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

all_peaks <- read.table("MSDIAL/Processing/Height_0_2019924102.txt", skip = 4, 
                        sep = "\t", header = T)[,c(1, 2, 3, 4, 29, 30, 31, 32)]

boxplot(log10(as.matrix(all_peaks[,c(5,6,7,8)])))

known_peaks <- filter(all_peaks, Metabolite.name!="Unknown")


mzml_files <- list.files("mzMLs", full.names = T)[c(1,5,6,7)]
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
