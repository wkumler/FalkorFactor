# Collects structural best-guesses from CSI:FingerID and compares them with
# metfRag and other software?

library(tidyverse)
library(rvest)

getSynonymsName <- function(cmpd_name, how_many=20){
  prolog <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
  tail <- "/cids/TXT"
  url <- paste0(prolog, cmpd_name, tail)
  CID <- readLines(url)
  getSynonymsCID(CID, how_many)
}

getSynonymsCID <- function(CID, how_many=20){
  prolog <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
  tail <- "/synonyms/TXT"
  url <- paste0(prolog, CID, tail)
  content <- head(readLines(url), how_many)
  paste0(content, collapse = "; ")
}

cmpd_file <- "XCMS/data_intermediate/sirius_project/output_dir/compound_identifications_adducts.tsv"
cmpd_raw <- read_tsv(cmpd_file, col_types = 
         cols(
           rank = col_double(),
           `#adducts` = col_double(),
           `#predictedFPs` = col_double(),
           Confidence_Score = col_character(),
           `CSI:FingerID_Score` = col_double(),
           Zodiac_Score = col_double(),
           TreeIsotope_Score = col_double(),
           molecularFormula = col_character(),
           adduct = col_character(),
           InChIkey2D = col_character(),
           InChI = col_character(),
           name = col_character(),
           smiles = col_character(),
           xlogp = col_logical(),
           pubchemids = col_character(),
           links = col_character(),
           dbflags = col_double(),
           id = col_character()
         ))
cmpd_data <- cmpd_raw %>%
  mutate(feature=sub(pattern = "[[:digit:]]*_", replacement = "", x = .$id)) %>%
  mutate(feature=sub(pattern = "_FEATURE1", replacement = "", x = .$feature)) %>%
  #mutate(top_cid=sub(pattern = ";.*", replacement = "", x=.$pubchemids)) %>%
  head() %>%
  #mutate(synonyms=sapply(X = .$top_cid, getSynonymsCID))
  mutate(synonyms=sapply(X = .$name, getSynonymsName))
cmpd_data %>% select(name, synonyms)
