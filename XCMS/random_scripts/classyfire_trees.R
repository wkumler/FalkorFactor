library(data.tree)
library(tidyverse)
library(tidygraph)


# Build tree ----
# Download raw data in JSON format
classyclasses <- jsonlite::fromJSON("http://classyfire.wishartlab.com/tax_nodes.json")
# Drop useless columns
classyclasses <- classyclasses[-1,c(2,3,1)]
# Rename to avoid data.tree naming problems
names(classyclasses) <- c("id", "parent", "class")
head(classyclasses)
# Save to disk
write.csv(classyclasses, "classyclasses.csv")

sirius_output <- "XCMS/data_intermediate/sirius_project/output_dir/" %>%
  paste0("canopus_summary.tsv") %>%
  readLines() %>%
  strsplit(split = "\t")
header <- sirius_output[[1]]
sirius_output <- sirius_output[-1] %>%
  do.call(what = rbind) %>%
  as.data.frame() %>%
  `names<-`(header)

classes_found <- sirius_output$`all classifications` %>%
  strsplit(split = "; ") %>%
  unlist() %>%
  table() %>%
  sort() %>%
  data.frame() %>%
  `names<-`(c("class", "count"))

compound_count <- merge(classyclasses, classes_found, by="class")
classyfire_tree <- FromDataFrameNetwork(classyclasses)
print(classyfire_tree, "class", "count")



