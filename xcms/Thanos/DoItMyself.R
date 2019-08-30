# Setup things ----
library(xcms)
library(tidyverse)

load("xcms/raw_data")
x <- raw_data %>%
  filterMsLevel(msLevel. = 1L) %>%
  selectFeatureData(fcol = c(MSnbase:::.MSnExpReqFvarLabels, "centroided")) %>%
  lapply(1:length(fileNames(.)), FUN=filterFile, object = .) %>%
  `[[`(1) %>%
  spectra()


# Generate peaklist ----
mzs <- lapply(x, mz)
mz <- unlist(mzs, use.names = FALSE)
int <- unlist(lapply(x, intensity), use.names = FALSE)
rts <- unlist(lapply(x, rtime))
rt <- rep(rts, sapply(mzs, length))
all_data <- data.frame(mz, int, rt)
ppm <- 2.5

# Optionally, remove all orphan data points
seq_data <- subset(all_data, duplicated(mz) | duplicated(mz, fromLast=TRUE))

# Determine typical m/z spread for a couple good peaks
# Unnecessary because m/z will be determined on-the-fly by ROI calculation, but looks nice
lmaoPlotEm <- function(data_i, default_layout=T) {
  Da_spread <- data_i$mz[which.max(data_i$int)]*ppm/1000000
  if(default_layout){
    layout(matrix(c(1,2), nrow = 2))
  }
  par(mar=c(0.1, 4.1, 2.1, 0.1))
  int_colors <- hcl.colors(100, palette = "plasma")[cut(data_i$int, breaks = 100)]
  plot(data_i$rt, data_i$mz, col=int_colors, xaxt="n", xlab="", pch=19, cex=1,
       ylim=c(data_i$mz[which.max(data_i$int)]-round(Da_spread, 5), 
              data_i$mz[which.max(data_i$int)]+round(Da_spread, 5)))
  legend("topleft", legend = paste("Min m/z:", round(min(data_i$mz), 5)))
  legend("topright", legend = paste("Max m/z:", round(max(data_i$mz), 5)))
  legend("bottomleft", legend = paste("m/z diff:", round(max(data_i$mz)-min(data_i$mz), 5)))
  legend("bottomright", legend = paste("epsilon:", round(Da_spread*2, 5)))
  par(mar=c(4.1, 4.1, 0.1, 0.1))
  plot(data_i$rt, data_i$int, col=int_colors, pch=19)
  if(default_layout){
    layout(1)
  }
}
seq_data %>% filter(mz>90.0910&mz<90.0925) %>% lmaoPlotEm()
seq_data %>% filter(mz>100.022&mz<100.025) %>% lmaoPlotEm() # WHAT IS HAPPENING WITH THIS ONE
seq_data %>% filter(mz>285.04&mz<285.06) %>% lmaoPlotEm()
seq_data %>% filter(mz>800.810&mz<800.815) %>% lmaoPlotEm()

# Generate ROI list via adapted ADAP algorithm
# Sort by intensity
ppm <- 2.5
peakwidth <- c(20, 80)
prefilter <- c(3, 100)

roi_list <- list()

if(exists("seq_data")){data <- seq_data}else{data <- all_data}
data <- filter(data, mz<70)

while(nrow(data)>0){
  point_of_interest <- data[which.max(data$int),]
  epsilon_Da <- point_of_interest$mz*ppm/1000000
  upper_roi_mz <- point_of_interest$mz+epsilon_Da
  lower_roi_mz <- point_of_interest$mz-epsilon_Da
  roi <- filter(data, mz>lower_roi_mz & mz<upper_roi_mz)
  
  # If the ROI can't contain a peak bc too short, remove it
  if(nrow(roi)<peakwidth[1]){
    data <- filter(data, mz<lower_roi_mz | mz>upper_roi_mz)
    next
  }
  
  # If there aren't enough points above the prefilter threshold
  runs_above_prefilter <- rle(roi$int>prefilter[2])
  max_run_length <- max(runs_above_prefilter$lengths[runs_above_prefilter$values])
  if(max_run_length<prefilter[2]){
    # Remove the ROI
    data <- filter(data, mz<lower_roi_mz | mz>upper_roi_mz)
    v <- roi
    next
  }
  
  roi_list[[length(roi_list)+1]] <- roi
  
  data <- filter(data, mz<lower_roi_mz | mz>upper_roi_mz)
}
