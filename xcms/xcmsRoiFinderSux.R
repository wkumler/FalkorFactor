# How does xcms find regions of interest???
# Does it really drop the ROI after a single missed scan
# Does it depend on intensity at all?
masses_of_interest <- c(90.0910, 90.0925)
# masses_of_interest <- c(117.0545, 117.0550) #Good ex of many ROIs and ragged EIC
# masses_of_interest <- c(100.023, 100.025)
# masses_of_interest <- c(117.048, 117.049)
# masses_of_interest <- c(285.040, 285.050)
# masses_of_interest <- c(800.810, 800.815) # Good ex findROI FUBAR

prefilter = c(1,1)
# prefilter = c(3, 10000) #Use this and below with above ex to show prefilter[1] works
# prefilter = c(10, 10000)

xlim <- NULL


scanindex <- as.integer(c(0, val_count[-length(val_count)]))
scantime = rt
scanrange <- c(1, length(scantime))
mz_span <- c(0.0005) # Maximum spread of m/z values across a well-defined peak, plus some buffer
ppm <- ceiling((mz_span*1000000)/132)
peakwidth <- c(20, 80)
min_peak_width <- min(peakwidth)/2
#min_centroids <- max(4, min_peak_width - 2)
min_centroids <- 9 #To parallel milliWave
noise = 0

roi_list <- .Call("findmzROI", 
                  mz, int, scanindex, 
                  as.double(c(0, 0)), 
                  as.integer(scanrange), 
                  as.integer(length(scantime)), 
                  as.double(ppm * 1e-06), 
                  as.integer(min_centroids), 
                  as.integer(prefilter), 
                  as.integer(noise), PACKAGE = "xcms")


xcms_found_rois <- as.data.frame(do.call(rbind, roi_list))
xcms_given_roi <- xcms_found_rois %>% 
  filter(mzmin>masses_of_interest[1]&mzmax<masses_of_interest[2])
will_same_roi <- all_data %>% 
  filter(mz>masses_of_interest[1]&mz<masses_of_interest[2])
runthings <- rle(rt%in%will_same_roi$rt)



layout(matrix(c(1,2,2,2), nrow = 4))
int_colors <- rev(hcl.colors(100, palette = "Terrain"))[cut(will_same_roi$int, breaks = 100)]
par(mar=c(0, 4.1, 0.1, 0.1))
plot(will_same_roi$rt, will_same_roi$int, pch=19, cex=1, col=int_colors, 
     xlim = xlim, ylab="Intensity", xlab="", xaxt="n")
par(mar=c(4.1, 4.1, 0, 0.1))
plot(will_same_roi$rt, will_same_roi$mz, pch=19, cex=1, col=int_colors,
     ylim=c(min(will_same_roi$mz)*0.999998, max(will_same_roi$mz)*1.000001),
     xlim=xlim, xlab = "Retention time (s)", ylab= "m/z")
if(nrow(xcms_given_roi)){
  for(row in 1:nrow(xcms_given_roi)){
    rect(xleft = rt[unlist(xcms_given_roi[row, "scmin"])],
         xright = rt[unlist(xcms_given_roi[row, "scmax"])],
         ytop = xcms_given_roi[row, "mzmax"],
         ybottom = xcms_given_roi[row, "mzmin"],
         border = "red", lwd=2)
  }
}
for(i in 1:length(runthings$lengths[c(T, F)])){
  if(runthings$values[1]){
    x0 <- rt[sum(runthings$lengths[0:(i*2-2)])]
    if(!length(x0)){
      x0 <- 1
    }
    x1 <- rt[sum(runthings$lengths[0:(i*2-1)])]
  } else {
    x0 <- rt[sum(runthings$lengths[1:(i*2-1)])]
    if(!length(x0)){
      x0 <- 1
    }
    x1 <- rt[sum(runthings$lengths[1:(i*2)])]
  }
  arrows(x0 = x0, x1 = x1,
         y0 = min(will_same_roi$mz)*0.999999,
         y1 = min(will_same_roi$mz)*0.999999,
         code = 3, angle = 90)
}

