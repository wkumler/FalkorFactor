# Peakfinding with different noise estimation algorithms

library(xcms)
library(tidyverse)

peak_height <- 10000
peak_width <- 50
scan_length <- 500
noise_sd <- 100
noise_mean <- 1000
noise_window_width <- 1 # As a multiple of peak_width


peak_values <- c(rep(0, (scan_length-peak_width)/2), 
                 table(cut(rnorm(peak_height/4*peak_width), breaks = peak_width)), 
                 rep(0, (scan_length-peak_width)/2)) %>%
  `+`(rnorm(n = scan_length, sd = noise_sd, mean = noise_mean))

max_index <- which.max(peak_values)
peak_start_index <- max_index-peak_width/2
peak_end_index <- max_index+peak_width/2
noise_window_indices <- list(c(peak_start_index-noise_window_width*peak_width, 
                               peak_start_index),
                             c(peak_end_index, 
                               peak_end_index+noise_window_width*peak_width))

left_noise_window_indices <- noise_window_indices[[1]][1]:noise_window_indices[[1]][2]
left_noise_window_indices <- left_noise_window_indices[left_noise_window_indices>0]
right_noise_window_indices <- noise_window_indices[[2]][1]:noise_window_indices[[2]][2]
right_noise_window_indices <- right_noise_window_indices[right_noise_window_indices<scan_length]


sd_noise_estimate <- sd(c(peak_values[left_noise_window_indices],
                          peak_values[right_noise_window_indices]),
                        na.rm = T)
median_noise_estimate <- median(c(peak_values[left_noise_window_indices],
                                  peak_values[right_noise_window_indices]),
                                na.rm = T)

plot(peak_values, pch=19)

abline(h=max(peak_values), col = "blue")
mtext(text = "Max\npeakheight", side = 4, at = max(peak_values), 
      col = "blue", line = 1)
arrows(x0 = noise_window_indices[[1]][1], x1 = noise_window_indices[[1]][2],
       y0 = median_noise_estimate, y1 = median_noise_estimate,
       col = "blue", angle = 90, code = 3, lwd = 2)
arrows(x0 = noise_window_indices[[2]][1], x1 = noise_window_indices[[2]][2],
       y0 = median_noise_estimate, y1 = median_noise_estimate,
       col = "blue", angle = 90, code = 3, lwd = 2)
text("Noise window", x = noise_window_indices[[2]][2], 
     y = median_noise_estimate+peak_height/10, col = "blue")

legend("topleft", legend = c("Signal-to-noise ratio",
                             paste("Max/Median:", round(max(peak_values)/median_noise_estimate)), 
                             paste("Max/SD:", round(max(peak_values)/sd_noise_estimate))),
       text.col = c("black", "green", "red"))
