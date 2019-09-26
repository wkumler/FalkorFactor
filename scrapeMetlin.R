# Getting metabolites database into reasonable shape to try peak identification


# Setup things ----
library(httr)
library(xml2)
library(dplyr)
pm <- function(x, d){c(x-d, x+d)}
mzr <- function(mz, ppm=2.5){round(pm(mz, mz*ppm/1000000), digits = 5)}

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
getMetlin <- function(cmpd_mz){
  mz_range <- mzr(cmpd_mz)
  metlin_url <- paste0("https://metlin.scripps.edu/advanced_search_result.php?",
                       "molid=&mass_min=", min(mz_range),
                       "&mass_max=", max(mz_range), "&Amino",
                       "Acid=add&drug=add&toxinEPA=add&keggIDFilter=add")
  metlin_raw <- GET(metlin_url)
  stop_for_status(metlin_raw)
  
  metlin_data <- metlin_raw %>%
    content() %>%
    xml_find_all(xpath = "//tbody") %>%
    xml_find_all(xpath="//td") %>%
    xml_text() %>%
    matrix(ncol=7, byrow=T) %>%
    as.data.frame() %>%
    `names<-`(c("exact_mass", "cmpd_name", "formula", 
                "CAS", "KEGG", "MSMS", "Structure"))
  
  print(paste("Metlin returned", nrow(metlin_data), "compounds between", 
        min(mz_range), "and", max(mz_range), "m/z"))
  print(paste("Of those,", sum(metlin_data$MSMS!="NO"), "have MS/MS data"))
  
  return(metlin_data)
}

sample_data <- getMetlin(135.054495)

metlin_data <- lapply(amino_masses, getMetlin)




test_url <- 'https://httpbin.org/user-agent'
metlin_url <- paste0("https://metlin.scripps.edu/advanced_search_result.php?",
                     "molid=&mass_min=", min(mz_range),
                     "&mass_max=", max(mz_range), "&Amino",
                     "Acid=add&drug=add&toxinEPA=add&keggIDFilter=add")
GET(url = test_url, add_headers('User-Agent'=user_agent))
GET(url = metlin_url, add_headers("User-Agent"=user_agent))
