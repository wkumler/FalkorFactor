

# Setup things ----
library(xcms)
library(dplyr)
library(pbapply)
setMSnbaseVerbose(TRUE)
register(SerialParam())



# File extraction and MSnExp creation ----
pooled_files <- list.files("mzMLs", pattern = "Full[1-3].mzML", full.names = TRUE)
sample_files <- list.files("mzMLs", pattern = "Smp", full.names = TRUE)
sample_meta <- gsub(x = basename(sample_files), replacement = "", 
                    pattern = "190715_Smp_FK180310") %>%
  data.frame(sample_name = .,
             sample_group = rep(c(rep("25m", 3), rep("DCM", 3)), 4),
             spindir = rep(c("Cyclonic", "Anti"), each=12),
             stringsAsFactors = FALSE) %>%
  AnnotatedDataFrame()
# raw_data <- readMSData(files = sample_files, msLevel. = 1, verbose = TRUE, 
#                        centroided. = TRUE, smoothed. = FALSE, pdata = sample_meta)
# save(raw_data, file = "XCMS/raw_data_MSnExp")


# Initial overview ----
load("XCMS/raw_data_MSnExp")
raw_data_short <- filterRt(raw_data, rt=c(120, 1100))
base_chrom <- chromatogram(raw_data_short, aggregationFun = "max")
save(base_chrom, file = "XCMS/base_chrom")
load("XCMS/base_chrom")
# DCM green, 25m blue
plot(base_chrom, col=c("#0000FF33", "#00FF0033")[as.factor(raw_data$sample_group)])
# Cyclonic red, Anti blue
plot(base_chrom, col=c("#FF000033", "#0000FF33")[as.factor(raw_data$spindir)])

library(pheatmap)
bpis_bin <- bin(base_chrom, binSize = 2)
cormat <- cor(log2(do.call(cbind, lapply(bpis_bin, intensity))))
colnames(cormat) <- rownames(cormat) <- raw_data$sample_name
ann <- data.frame(group = raw_data$sample_group)
rownames(ann) <- raw_data$sample_name
pheatmap(cormat, annotation = ann)

