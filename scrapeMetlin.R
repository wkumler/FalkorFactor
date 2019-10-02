# Obtain data from Metlin by submitting GET request behind their Javascript wall
# Three main functions: getMetlinMz, getMetlinName, getMetlinMS2


# Setup things ----
library(httr)
library(xml2)
library(dplyr)
pm <- function(x, d){c(x-d, x+d)}
mzr <- function(mz, ppm){round(pm(mz, mz*ppm/1000000), digits = 5)}

# getMetlin functions ----
getMetlinMz <- function(cmpd_mz, ppm=2.5){
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
  
  MSMS_subset <- subset(search_data, MSMS=="experimental")
  
  print(paste("Metlin returned", nrow(search_data), "compound(s) between", 
        min(mz_range), "and", max(mz_range), "m/z with", 
        length(levels(search_data$formula)), "unique formula(s):", 
        paste(levels(search_data$formula), collapse = ", ")))
  if(nrow(MSMS_subset)){
    print(paste("Of those,", nrow(MSMS_subset), 
                "have experimental MS/MS data:", 
                paste(as.character(MSMS_subset$cmpd_name), collapse = ", ")))
  } else {
    print("None of these compounds have experimental MS/MS data")
  }
  
  return(search_data)
}

getMetlinName <- function(name){
  metlin_data <- paste0("https://metlin.scripps.edu/advanced_search_result.php?", 
                        "name=", name, "&AminoAcid=add&drug=add&toxinEPA=add&",
                        "keggIDFilter=add") %>%
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
  
  MSMS_subset <- subset(search_data, MSMS=="experimental")
  
  print(paste("Metlin returned", nrow(search_data), "compound(s) with name", 
              name, " with", 
              length(levels(search_data$formula)), "unique formula(s):", 
              paste(levels(search_data$formula), collapse = ", ")))
  if(nrow(MSMS_subset)){
    print(paste("Of those,", nrow(MSMS_subset), 
                "have experimental MS/MS data:", 
                paste(as.character(MSMS_subset$cmpd_name), collapse = ", ")))
  } else {
    print("None of these compounds have experimental MS/MS data")
  }
  
  return(search_data)
}

getMetlinMS2 <- function(cmpd_id){
  metlin_data <- paste0("https://metlin.scripps.edu/showChart.php?molid=", 
                        cmpd_id,"&etype=experimental") %>%
    GET() %>% stop_for_status()
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
  ms2_polarities <- ifelse(grepl(ms2_titles, pattern = "(\\+)"), "+", "-")
  ms2_voltages <- gregexpr("\\d+\\.?\\d*", ms2_titles) %>%
    regmatches(x=ms2_titles) %>%
    unlist() %>% as.numeric()
  ms2_adducts <- gregexpr("\\[M.*\\]", ms2_titles) %>%
    regmatches(x = ms2_titles) %>% 
    unlist()
  
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
  for(i in seq_along(ms2_list)){
    ms2_list[[i]] <- cbind(ms2_polarities[i], ms2_adducts[i], 
                           ms2_voltages[i], ms2_list[[i]])
  }
  ms2_df <- as.data.frame(do.call(rbind, ms2_list))
  names(ms2_df) <- c("polarity", "adduct", "voltage", "frag_mass", "frag_int")
  
  print(paste("Metlin had", length(ms2_titles), "MS2 records for this compound,",
              "with collision energies of", 
              paste(paste0(ms2_polarities, ms2_voltages), collapse=", ")))
  
  return(ms2_df)
}

# Using getMetlin functions ----

sample_data <- getMetlinMz(117.078979)
sample_data <- getMetlinMz(117.078979, ppm = 500)
sample_data <- getMetlinName("Arsenobetaine")

sample_ms2_cmpd <- sample_data %>% filter(MSMS=="experimental") %>% slice(1)

sample_ms2 <- getMetlinMS2(sample_ms2_cmpd$cmpd_id)


# Visualize the ill-gotten data ----
split_volt_pos_ms2 <- sample_ms2 %>%
  subset(polarity=="+") %>%
  split(.$voltage)
layout(matrix(c(1, rep(2:(length(split_volt_pos_ms2)+1), each=2),  1), ncol = 1))
par(mar=c(0.1, 4.1, 0.1, 0.1))
plot.new()
text(x = 0.5, y=1, labels = sample_ms2_cmpd$cmpd_name, cex=3)
for(i in split_volt_pos_ms2){
  plot(i$frag_mass, i$frag_int, xlab = "", 
       ylab=paste("Voltage", unique(i$voltage)),
       xlim=c(0, max(sample_ms2$frag_mass)), 
       xaxt="n", yaxt="n", type="n", ylim=c(0, 120))
  segments(x0 = i$frag_mass, x1 = i$frag_mass,
           y0 = 0, y1 = i$frag_int)
  axis(side = 2, at = c(0, 50, 100), labels = c(0, 50, 100))
}
axis(side = 1)

library(ggplot2)
gp <- sample_ms2 %>%
  filter(polarity=="+") %>%
    #filter(voltage==20) %>%
    ggplot(label=frag_mass) +
    geom_segment(aes(yend=0, x=frag_mass, y=frag_int, xend=frag_mass)) +
    #geom_hline(yintercept=0) +
    facet_wrap(~voltage, ncol = 1) +
    theme_bw() +
    xlim(0, max(sample_ms2$frag_mass))
gp

library(plotly)
ggplotly(gp, tooltip = c("y", "x"))


# Begin preparations to scrape the whole database ----
getMetlinMzRange <- function(mz_min, mz_max){
  metlin_data <- paste0("https://metlin.scripps.edu/advanced_search_result.php?",
                        "molid=&mass_min=", mz_min,
                        "&mass_max=", mz_max, "&Amino",
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
  return(search_data)
}

for(i in seq(171, 180)){
  print(c(i, i+0.06))
  print(nrow(getMetlinMzRange(i, (i+0.06))))
  print(c(i+0.06, i+0.1))
  print(nrow(getMetlinMzRange(i+0.06, (i+.1))))
  print(c(i+0.1, i+1))
  print(nrow(getMetlinMzRange(i+0.1, i+1)))
  Sys.sleep(runif(1))
}
dput(c(0.06, 0.1, 1)+rep(171:180, each=3))

c(0, 59, 72, 80, 86, 90, 95, 98, 100, 102, 104, 106, 108, 110, 111:124,
  124.1, 125, 125.1, 126, 126.1, 127, 127.1, 128, 128.08, 129, 129.08,
  130, 130.08, 131, 131.08, 132, 132.08, 133, 133.08, 134, 134.08, 
  135, 135.08, 136, 136.08, 137, 137.08, 138, 138.08, 139, 139.08, 
  140, 140.08, 140.12, 141, 141, 141.06, 141.1, 142, 142.06, 142.1, 143, 
  143.06, 143.1, 144, 144.06, 144.1, 145, 145.06, 145.1, 146,
  146.06, 146.1, 147, 147.06, 147.1, 148, 148.06, 148.1, 149, 
  149.06, 149.1, 150, 150.06, 150.1, 151,
  151.06, 151.1, 152, 152.06, 152.1, 153, 153.06, 153.1, 154, 
  154.06, 154.08, 154.1, 155, 155.06, 155.1, 156, 156.06, 156.1, 156.3, 157, 157.06, 
  157.1, 158, 158.06, 158.1, 159, 159.06, 159.1, 160, 160.06, 160.08, 160.1, 
  161, 161.06, 161.1, 162, 162.06, 162.1, 163, 163.06, 163.1, 164, 
  164.06, 164.08, 164.1, 165, 165.06, 165.1, 166, 166.06, 166.08, 166.1, 167, 167.06, 
  167.1, 168, 168.04, 168.06, 168.1, 168.3, 169, 169.06, 169.1, 170, 170.04,
  170.06, 170.1, 170.3, 171, 171.06, 171.1, 172, 172.06, 172.1, 172.3, 
  173, 173.06, 173.1, 174, 174.06, 174.1, 174.3, 175, 175.06, 175.08, 
  175.1, 176, 176.06, 176.08, 176.1, 177, 177.06, 177.1, 177.3, 178)
