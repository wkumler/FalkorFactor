# Getting metabolites database into reasonable shape to try peak identification

library(httr)
library(xml2)
pm <- function(x, d){c(x-d, x+d)}
mzr <- function(mz, ppm=2.5){round(pm(mz, mz*ppm/1000000), digits = 5)}

getMetlin <- function(cmpd_mz){
  mz_range <- mzr(cmpd_mz)
  metlin_url <- paste0("https://metlin.scripps.edu/advanced_search_result.php?",
                       "molid=&mass_min=", min(mz_range),
                       "&mass_max=", max(mz_range), "&Amino",
                       "Acid=add&drug=add&toxinEPA=add&keggIDFilter=add")
  metlin_data <- metlin_url %>%
    read_html() %>%
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

amino_masses <- "http://www.matrixscience.com/help/aa_help.html" %>%
  read_html() %>%
  xml_find_all(xpath="//td") %>%
  xml_text() %>%
  matrix(ncol=6, byrow=T) %>%
  gsub("\r\n    ", "", .)
amino_names <- amino_masses[,1]
amino_masses <- as.numeric(amino_masses[,4])
names(amino_masses) <- amino_names
amino_masses <- amino_masses[!is.na(amino_masses)]+18.010565

lapply(amino_masses, getMetlin)
