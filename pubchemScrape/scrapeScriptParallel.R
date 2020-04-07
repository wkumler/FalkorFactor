# PubChem scraper

setwd("C:\\Users\\Owner\\Desktop\\Will\\pubchemScrape")

library(httr)
library(xml2)
library(dplyr)
library(BiocParallel)

baseurl <- "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/"
index <- GET(url = baseurl)
filenames <- xml_find_all(content(index), xpath = "//a") %>%
  `[`(grep(pattern = "\\.gz\"", x = .)) %>%
  gsub(pattern = '<a href="', replacement = "") %>%
  gsub(pattern = '\\\">C.*', replacement = "")

register(BPPARAM = SnowParam(workers = 12, tasks = length(filenames)))
all_formulae <- bplapply(filenames, function(file_name, baseurl){
  url <- paste0(baseurl, file_name)
  download.file(url, destfile = file_name)
  file_data <- readLines(file_name)
  formula_lines <- grep(pattern = "PUBCHEM_MOLECULAR_FORMULA", 
                        file_data, fixed = TRUE)
  file.remove(file_name)
  return(unique(file_data[formula_lines+1]))
}, baseurl)
formulae_unique <- unique(unlist(all_formulae))
saveRDS(formulae_unique, file = "unique_formulae.rds")