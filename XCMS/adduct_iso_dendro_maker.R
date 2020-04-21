library(vegan)
library(dendextend)
xdata_filled <- readRDS("XCMS/temp_data/current_xdata_filled.rds")

feature_values <- featureValues(xdata_filled)
feature_defs <- featureDefinitions(xdata_filled)


dend <- feature_values %>%
  decostand(method = 'standardize', 1, na.rm = T) %>%
  dist() %>% hclust(method='average') %>%
  as.dendrogram()
pdf(file = "XCMS/temp_data/dendro.pdf", height = 100, width = 8)
par(mar=c(4.1, 4.1, 4.1, 16.1))
plot(dend, horiz=TRUE)
dev.off()

titledata <- as.data.frame(feature_defs[,c("mzmed", "rtmed")]) %>%
  mutate(feature=rownames(.)) %>%
  mutate(title=paste0(feature, ": m/z=", round(mzmed, 5), ", RT=", round(rtmed))) %>%
  pull(title)
labels(dend) <- titledata[order(labels(dend))]
pdf(file = "XCMS/temp_data/dendro_addinfo.pdf", height = 100, width=8)
par(mar=c(4.1, 4.1, 4.1, 16.1))
plot(dend, horiz=TRUE)
dev.off()
