# Script for running xcms on Falkor data files

# Setup things ----
library(xcms)
library(tidyverse)



# Data import ----
mzxml_path <- "Z:/1_QEdata/LTC/DATA/HILIC/190718_DepthProfiles_FK180310/"
mzxml_pos_files <- list.files(paste0(mzxml_path, "positive"), full.names = T)
useful_files <- mzxml_pos_files[c(1,5:7,11:13,14:37)]

pdata <- data.frame(sample_name = sub(basename(useful_files), pattern = ".mzXML",
                                      replacement = "", fixed = TRUE),
                    sample_group = c("Seawater_filter_blank",
                                     rep("Full_pooled", 3),
                                     rep("Half_pooled", 3),
                                     rep("Sample", length(14:37))))

raw_data <- readMSData(files = useful_files, 
                       pdata = new("NAnnotatedDataFrame", pdata),
                       mode = "onDisk")
save(raw_data)
load(raw_data)

group_colors <- c(rgb(0,0,0), rgb(1,0,0,0.5), rgb(0,1,0,0.5), rgb(0,0,1,0.1))
names(group_colors) <- c("Seawater_filter_blank", "Full_pooled", 
                         "Half_pooled", "Sample")



# Initial data inspection ----

bpis <- chromatogram(raw_data, aggregationFun = "max")
plot(bpis, col = group_colors[raw_data$sample_group])

tc <- split(tic(raw_data), f = fromFile(raw_data))
boxplot(tc, col = group_colors[raw_data$sample_group],
        ylab = "intensity", main = "Total ion current")