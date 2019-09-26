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
  metlin_data <- paste0("https://metlin.scripps.edu/advanced_search_result.php?",
                       "molid=&mass_min=", min(mz_range),
                       "&mass_max=", max(mz_range), "&Amino",
                       "Acid=add&drug=add&toxinEPA=add&keggIDFilter=add") %>%
    GET() %>%
    stop_for_status() %>%
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

getMetlinMS2 <- function(cmpd_id){
  metlin_data <- paste0("https://metlin.scripps.edu/showChart.php?molid=", cmpd_id,
                        "&etype=experimental") %>%
    GET() %>%
    stop_for_status()
  ms2_raw <- metlin_data %>%
    content(encoding = "UTF-8") %>%
    xml_find_all(xpath = '/html/head/script[5]') %>%
    xml_text() %>%
    gsub(pattern = ".*series: \\[", replacement = "") %>%
    gsub(pattern = " ]}]\n                }); \n\n            });\n        ",
         replacement = "") %>%
    gsub(pattern = "&nbsp;", replacement = "") %>%
    strsplit(split = "\\{name: ") %>%
    unlist() %>%
    `[`(-1) %>%
    lapply(strsplit, split=",data:\\[")
  
  if(!length(ms2_raw)){stop("No MS2 data found. Are you sure Metlin has this data?")}
  
  ms2_titles <- sapply(ms2_raw, `[[`, 1)[1,] #Not sure why this works but OK
  ms2_polarities <- ifelse(grepl(ms2_titles, pattern = "(\\+)"), "Pos", "Neg")
  ms2_voltages <- as.numeric(unlist(regmatches(ms2_titles, m = gregexpr("\\d+\\.?\\d*", ms2_titles))))
  ms2_adducts <- unlist(regmatches(ms2_titles, m = gregexpr("\\[M.*\\]", ms2_titles)))
  
  ms2_list <- sapply(ms2_raw, `[[`, 1)[2,] %>%
    gsub(pattern = " ]},", replacement = "") %>%
    lapply(function(x){
      wo1 <- substring(x, 2)
      wo2 <- substring(wo1, 1, nchar(wo1))
      ms2_list <- unlist(strsplit(wo2, "},\\{"))
      ms2_list <- lapply(ms2_list, function(x)
        as.numeric(unlist(regmatches(x, m = gregexpr("\\d+\\.?\\d*", x)))))
      ms2_df <- as.data.frame(do.call(rbind, ms2_list))
      names(ms2_df) <- c("mass", "rel_intensity")
      return(ms2_df)
    })
  v <- mapply(FUN= list, ms2_polarities, ms2_adducts, ms2_voltages)
}
