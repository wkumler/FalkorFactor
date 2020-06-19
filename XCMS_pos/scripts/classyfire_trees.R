library(data.tree)
library(tidyverse)
library(tidygraph)


# Build tree ----
raw_json <- jsonlite::fromJSON("http://classyfire.wishartlab.com/tax_nodes.json")
raw_json <- raw_json[-1,c(2,3,1)]
names(raw_json) <- c("id", "parent", "class")
head(raw_json)

sirius_output <- "XCMS_pos/data_intermediate/sirius_project/output_dir/" %>%
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

compound_count <- merge(raw_json, classes_found, by="class")
classyfire_tree <- FromDataFrameNetwork(raw_json)
print(classyfire_tree, "class", "count")



