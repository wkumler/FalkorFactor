# Script for running xcms on Falkor data files

# Setup things ----
library(xcms)
library(tidyverse)
library(pheatmap)
will_plotXIC <- function (x, main = "", col = NA, colramp = topo.colors, 
                          grid.color = "lightgrey", pch = 21, yl = NULL,
                          layout = matrix(1:2, ncol = 1), mn = NULL, 
                          spectrumnum = NULL, ...) {
  x <- suppressWarnings(as(x, "data.frame"))
  if(!is.null(spectrumnum)){
    x <- x[x$file%in%spectrumnum,]
  }
  bpi <- unlist(lapply(split(x$i, x$rt), max, na.rm = TRUE))
  brks <- lattice::do.breaks(range(x$i), nint = 256)
  par(mar = c(0, 4.5, 2, 1))
  plot(as.numeric(names(bpi)), bpi, xaxt = "n", col = col, 
       main = mn, bg = lattice:::level.colors(bpi, at = brks, col.regions = colramp), 
       xlab = "", pch = pch, ylab = "", las = 2, 
       ...)
  mtext(side = 4, line = 0, "Intensity", cex = par("cex.lab"))
  grid(col = grid.color)
  par(mar = c(3.5, 4.5, 0, 1))
  plot(x$rt, x$mz, main = "", pch = pch, col = col, xlab = "", 
       ylab = "", yaxt = "n", ylim = yl,
       bg = lattice:::level.colors(x$i, at = brks, col.regions = colramp), ...)
  axis(side = 2, las = 2)
  grid(col = grid.color)
  mtext(side = 1, line = 2.5, "Retention time", cex = par("cex.lab"))
  mtext(side = 4, line = 0, "m/z", cex = par("cex.lab"))
}


# Data import ----
mzml_path <- "mzMLs"
mzml_pos_files <- list.files(mzml_path, full.names = T)
useful_files <- mzml_pos_files[c(1,5:7,17:40)]

pdata <- data.frame(sample_name = sub(basename(useful_files), pattern = ".mzML",
                                      replacement = "", fixed = TRUE),
                    sample_group = factor(c("Seawater_filter_blank",
                                            rep("Full_pooled", 3),
                                            rep("Sample", length(17:40))),
                                          levels = c("Seawater_filter_blank", "Full_pooled", "Sample")),
                    env_group = factor(c("Seawater_filter_blank",
                                         rep("Full_pooled", 3),
                                         rep(c(rep("25m", 3), rep("DCM", 3)), length(17:40)/6)),
                                       levels = c("Seawater_filter_blank", "Full_pooled", "25m", "DCM")),
                    loc_group = factor(c("Seawater_filter_blank",
                                  rep("Full_pooled", 3),
                                  c(rep("CCW", length(17:40)/2), rep("Clockwise", length(17:40)/2))),
                                  levels = c("Seawater_filter_blank", "Full_pooled", "CCW", "Clockwise")))

# raw_data <- readMSData(files = useful_files,
#                        pdata = new("NAnnotatedDataFrame", pdata),
#                        mode = "onDisk")
# save(raw_data, file = "xcms/raw_data")
load("xcms/raw_data")

sample_group_colors <- c(rgb(0,0,0), rgb(1,0,0,0.5), rgb(0,0,1,0.1)) %>%
  `names<-`(c("Seawater_filter_blank", "Full_pooled", "Sample"))
env_group_colors <- c(rgb(0,0,0), rgb(1,0,0,0.5), rgb(0,1,0,0.2), rgb(0,0,1,0.2)) %>%
  `names<-`(c("Seawater_filter_blank", "Full_pooled", "25m", "DCM"))
loc_group_colors <- c(rgb(0,0,0), rgb(1,0,0,0.5), rgb(0,1,0,0.2), rgb(0,0,1,0.2)) %>%
  `names<-`(c("Seawater_filter_blank", "Full_pooled", "Clockwise", "CCW"))


# Initial data inspection ----

bpis <- chromatogram(raw_data, aggregationFun = "max")
save(bpis, file = "xcms/bpis")
load("xcms/bpis")
plot(bpis, col = sample_group_colors[raw_data$sample_group])
plot(bpis, col = env_group_colors[raw_data$env_group])
legend("top", legend = c("Seawater_filter_blank", "Full_pooled", "25m", "DCM"),
       col = c(rgb(0,0,0), rgb(1,0,0), rgb(0,1,0), rgb(0,0,1)), lty = 1)
plot(bpis, col = loc_group_colors[raw_data$loc_group])
legend("top", legend = c("Seawater_filter_blank", "Full_pooled", "Clockwise", "CCW"),
       col = c(rgb(0,0,0), rgb(1,0,0), rgb(0,1,0), rgb(0,0,1)), lty = 1)



tc <- split(tic(raw_data), f = fromFile(raw_data))
layout(matrix(c(1,1,2,2,2,2), byrow = T))
par(mar=c(2.1, 4.1, 0.1, 0.1))
boxplot(tc, col = sample_group_colors[raw_data$sample_group],
        ylab = "intensity")
boxplot(tc, col = sample_group_colors[raw_data$sample_group],
        ylab = "intensity", ylim = c(0, 350000000))
layout(1)




## Calculate correlation on the log2 transformed base peak intensities
bpis_bin <- bin(bpis, binSize = 2)
cormat <- cor(log2(do.call(cbind, lapply(bpis_bin, intensity))))
colnames(cormat) <- rownames(cormat) <- raw_data$sample_group
pheatmap(cormat)



# Peakfinding using CentWave ----
stds <- read.csv("xcms/Ingalls_Lab_Standards_Will.csv", stringsAsFactors = F)

plotSTD <- function(std_row){
  par(mar=c(2.1, 2.1, 1.1, 0.1))
  rtr<-c(stds$rt.sec[std_row]-100, stds$rt.sec[std_row]+100)
  mzr<-c(stds$m.z[std_row]-0.01, stds$m.z[std_row]+0.01)
  chr_raw <- chromatogram(raw_data, mz = mzr, rt = rtr)
  plot(chr_raw, col = sample_group_colors[chr_raw$sample_group], main=stds$Compound.Name[std_row])
}

good_peak <- 18 # Found by cycling through standards one at a time

random_stds <- c(good_peak, sample(1:dim(stds)[1], 3))
par(mfrow=c(2,2))
for(i in random_stds){plotSTD(i)}
peakwidth <- c(20, 80) # Minimum and maximum width of a few random peaks, in seconds

par(mfrow=c(2,1))
raw_data %>%
  filterRt(rt = c(stds$rt.sec[good_peak]-100, stds$rt.sec[good_peak]+100)) %>%
  filterMz(mz = c(stds$m.z[good_peak]-0.001, stds$m.z[good_peak]+0.001)) %>%
  will_plotXIC(spectrumnum = 1:28)
par(mfrow=c(1,1))
mz_span <- c(0.0005) # Maximum spread of m/z values across a well-defined peak, plus some buffer
ppm <- ceiling((mz_span*1000000)/132)



# Define our CentWave parameters based on the above
cwp <- CentWaveParam(peakwidth = peakwidth, ppm = ppm, 
                     snthresh = 1, prefilter = c(3,10000))


# Plot a chromatogram for the given peak
good_chr_raw <- chromatogram(raw_data,
                        mz = c(stds$m.z[good_peak]-mz_span, 
                               stds$m.z[good_peak]+mz_span),
                        rt = c(stds$rt.sec[good_peak]-max(peakwidth), 
                               stds$rt.sec[good_peak]+max(peakwidth)))

# Check peakfinding with the new parameters
xchr <- findChromPeaks(good_chr_raw, param = cwp)
chromPeaks(xchr)
chromPeakData(xchr)
plot(xchr, peakType="rectangle", peakCol=sample_group_colors[xchr$sample_group])



# Find peaks in the full data set
# xdata <- findChromPeaks(raw_data, param = cwp) #Takes about 25 mins in serial on laptop
# save(xdata, file = "xcms/xdata")
load("xcms/xdata")

all_peaks <- chromPeaks(xdata)
head(chromPeaks(xdata))
tail(chromPeaks(xdata))

plotChromPeaks(xdata, file = 28)
par(mar=c(4.1, 10.1, 0.1, 0.1))
plotChromPeakImage(xdata)
par(mar=c(4.1, 4.1, 2.1, 0.1))
split(log2(chromPeaks(xdata)[, "into"]),
              f = chromPeaks(xdata)[, "sample"]) %>%
  boxplot(varwidth = TRUE, col=sample_group_colors[chr_ex$sample_group],
          ylab = expression(log[2]~intensity), 
          main = "Peak intensities")
grid(nx = NA, ny = NULL)



chr_ex <- chromatogram(xdata, 
                       mz = c(stds$m.z[good_peak]-0.001, stds$m.z[good_peak]+0.001),
                       rt = c(stds$rt.sec[good_peak]-100, stds$rt.sec[good_peak]+100))
chromPeaks(chr_ex)
plot(chr_ex, col = "black", 
     peakType = "rectangle", peakCol = sample_group_colors[chr_ex$sample_group],
     peakBg = NA)

chromPeaks(xdata, mz = c(stds$m.z[good_peak]-0.0005, stds$m.z[good_peak]+0.0005),
           rt = c(stds$rt.sec[good_peak]-100, stds$rt.sec[good_peak]+100))

chr_ex <- chromatogram(xdata, mz = c(stds$m.z[good_peak]-0.001, stds$m.z[good_peak]+0.001))
plot(chr_ex, peakType = "rectangle")
