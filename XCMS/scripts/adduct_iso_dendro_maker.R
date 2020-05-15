library(tidyverse)
library(vegan)
library(dendextend)
xdata_filled <- readRDS("XCMS/temp_data/current_xdata_filled.rds")

feature_values <- xcms::featureValues(xdata_filled) %>%
  #`[`(,grepl(pattern = "Smp", x = colnames(.)))
  `[`(,-c(1:7, 32:39, 41))
feature_defs <- xcms::featureDefinitions(xdata_filled)


dend <- feature_values %>%
  decostand(method = 'standardize', 1, na.rm = T) %>%
  dist() %>% hclust(method='average') %>%
  as.dendrogram()

vegdend <- feature_values %>%
  decostand(method = 'standardize', MARGIN = 1, na.rm = TRUE) %>%
  vegdist(method = "euclidian", na.rm = TRUE) %>% hclust(method='average') %>%
  as.dendrogram()


titledata <- as.data.frame(feature_defs[,c("mzmed", "rtmed")]) %>%
  mutate(feature=rownames(.)) %>%
  mutate(title=paste0(feature, ": m/z=", round(mzmed, 5), ", RT=", round(rtmed))) %>%
  pull(title)
labels(vegdend) <- titledata[order(labels(vegdend))]
pdf(file = "XCMS/temp_data/dendro_addinfo.pdf", height = 100, width=8)
par(mar=c(4.1, 4.1, 4.1, 16.1))
plot(vegdend, horiz=TRUE)
dev.off()



