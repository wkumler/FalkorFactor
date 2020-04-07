# PubChem scraper

setwd("G://My Drive//FalkorFactor/pubchemScrape")

library(httr)
library(xml2)
library(dplyr)

baseurl <- "https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/CURRENT-Full/SDF/"
index <- GET(url = baseurl)
filenames <- xml_find_all(content(index), xpath = "//a") %>%
  `[`(-1) %>%
  `[`(-length(.)) %>%
  `[`(-grep(pattern = ".md5", x = .)) %>%
  gsub(pattern = '<a href="', replacement = "") %>%
  gsub(pattern = '\\\">C.*', replacement = "")

all_formulae <- list()
pb <- txtProgressBar(min = 0, max = length(filenames))
for(file_name in filenames){
  url <- paste0(baseurl, file_name)
  download.file(url, destfile = "current_pubchem_data.gz")
  file_data <- readLines("current_pubchem_data.gz")
  formula_lines <- grep(pattern = "PUBCHEM_MOLECULAR_FORMULA", 
                        file_data, fixed = TRUE)
  formulae <- unique(file_data[formula_lines+1])
  setTxtProgressBar(pb, which(filenames==file_name))
}
formulae_unique <- unique(unlist(all_formulae))
saveRDS(formulae_unique, file = "unique_formulae.rds")