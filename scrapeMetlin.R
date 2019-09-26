# Getting metabolites database into reasonable shape to try peak identification


# Setup things ----
library(httr)
library(xml2)
library(dplyr)
pm <- function(x, d){c(x-d, x+d)}
mzr <- function(mz, ppm){round(pm(mz, mz*ppm/1000000), digits = 5)}

amino_masses <- "http://www.matrixscience.com/help/aa_help.html" %>%
  GET() %>%
  content() %>%
  xml_find_all(xpath="//td") %>%
  xml_text() %>%
  matrix(ncol=6, byrow=T) %>%
  gsub("\r\n    ", "", .) %>%
  `rownames<-`(.[,1]) %>%
  `[`(1:24, 4) %>%
  `mode<-`("numeric") %>% 
  na.omit() %>%
  `+`(18.010565)

# getMetlin code ----
getMetlin <- function(cmpd_mz, ppm=2.5){
  mz_range <- mzr(cmpd_mz, ppm)
  metlin_url <- paste0("https://metlin.scripps.edu/advanced_search_result.php?",
                       "molid=&mass_min=", min(mz_range),
                       "&mass_max=", max(mz_range), "&Amino",
                       "Acid=add&drug=add&toxinEPA=add&keggIDFilter=add")
  metlin_raw <- GET(metlin_url)
  stop_for_status(metlin_raw)
  
  metlin_data <- metlin_raw %>%
    content(encoding = "UTF-8") %>%
    xml_find_all(xpath = "//tbody")
  
  search_ids <- metlin_data %>%
    xml_find_all(xpath = "//th[@scope]/a") %>%
    xml_text()
  
  search_data <- metlin_data %>%
    xml_find_all(xpath="//td") %>%
    xml_text() %>%
    matrix(ncol=7, byrow=T) %>%
    cbind(search_ids, .) %>%
    as.data.frame() %>%
    `names<-`(c("cmpd_id", "exact_mass", "cmpd_name", "formula", 
                "CAS", "KEGG", "MSMS", "Structure"))
  
  print(paste("Metlin returned", nrow(search_data), "compound(s) between", 
        min(mz_range), "and", max(mz_range), "m/z with", 
        length(levels(search_data$formula)), "unique formula(s):", 
        paste(levels(search_data$formula), sep = ", ")))
  print(paste("Of those,", sum(search_data$MSMS!="NO"), "have MS/MS data"))
  
  return(search_data)
}

sample_data <- getMetlin(135.054495)
sample_data <- getMetlin(117.078979)

metlin_data <- lapply(amino_masses[1:3], getMetlin)



# getHMDB code ----

