# Script to run SIRIUS (4.4.17), Rdisop, and check isotope signatures and obtain
# high-confidence formula assignments to features
# Expects a feature list (usually created by peakpicking.R)
# Isotope checking expects a peak list (usually created by peakpicking.R)
# Needs bugtesting after WD restructuring


#Only for RScript running
setwd(r"(G:\My Drive\FalkorFactor)")



# Setup things ----
library(tidyverse)
library(data.table)
library(pbapply)
# Requires mzR to load MSMS data
final_features <- read.csv(file = "XCMS/data_pretty/final_features.csv")
final_peaks <- read.csv(file = "XCMS/data_pretty/final_peaks.csv")
sirius_project_dir <- "XCMS/data_intermediate/sirius_project"




# Functions ----
pmppm <- function(mass, ppm=4){c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))}
formula2elements <- function(formula_vec){
  split_formulas <- formula_vec %>%
    gregexpr(pattern = "[A-Z][a-z]*[0-9]*") %>%
    regmatches(x = formula_vec)
  elements_in <- lapply(split_formulas, gsub, pattern = "[0-9]", replacement = "")
  element_counts <- lapply(split_formulas, gsub, pattern = "[A-z]", replacement = "")
  element_counts <- lapply(element_counts, function(x){
    x[!nchar(x)]<-1
    return(as.numeric(x))
  })
  mapply(`names<-`, element_counts, elements_in, SIMPLIFY = FALSE)
}
grabSingleFileMS2 <- function(filename){
  msdata <- mzR::openMSfile(filename)
  fullhd <- mzR::header(msdata)
  ms2rows <- seq_len(nrow(fullhd))[fullhd$msLevel>1]
  spectra_list <- lapply(ms2rows, function(x){
    rtime <- fullhd[x, "retentionTime"]
    premz <- fullhd[x, "precursorMZ"]
    fragments <- mzR::peaks(msdata, x)
    return(cbind(rtime, premz, fragments))
  })
  all_data <- `names<-`(as.data.frame(do.call(rbind, spectra_list)), 
                        c("rt", "premz", "fragmz", "int"))
  return(all_data)
}
mgf_maker <- function(feature_msdata, ms1, ms2, output_file){
  if(!nrow(ms2)){
    outtext <- c("BEGIN IONS",
                 paste0("PEPMASS=", feature_msdata$mzmed),
                 "MSLEVEL=1",
                 "CHARGE=1+",
                 apply(ms1, 1, paste, collapse=" "),
                 "END IONS",
                 "")
  } else {
    outtext <- c("BEGIN IONS",
                 paste0("PEPMASS=", feature_msdata$mzmed),
                 "MSLEVEL=1",
                 "CHARGE=1+",
                 apply(ms1, 1, paste, collapse=" "),
                 "END IONS",
                 "",
                 "BEGIN IONS",
                 paste0("PEPMASS=", feature_msdata$mzmed),
                 "MSLEVEL=2",
                 "CHARGE=1+",
                 apply(ms2[ms2$voltage==20, c("fragmz", "int")], 
                       1, paste, collapse=" "),
                 "END IONS",
                 "",
                 "BEGIN IONS",
                 paste0("PEPMASS=", feature_msdata$mzmed),
                 "MSLEVEL=2",
                 "CHARGE=1+",
                 apply(ms2[ms2$voltage==35, c("fragmz", "int")], 
                       1, paste, collapse=" "),
                 "END IONS",
                 "BEGIN IONS",
                 paste0("PEPMASS=", feature_msdata$mzmed),
                 "MSLEVEL=2",
                 "CHARGE=1+",
                 apply(ms2[ms2$voltage==50, c("fragmz", "int")], 
                       1, paste, collapse=" "),
                 "END IONS")
  }
  writeLines(outtext, con = output_file)
}
isocheck <- function(feature_num, final_peaks=final_peaks, printplot=FALSE){
  ft_isodata <- final_peaks %>% filter(feature==feature_num) %>%
    select(M_area, C13_area, N15_area, O18_area, X2C13_area, S34_area, S33_area) %>%
    pivot_longer(cols = starts_with(c("C", "X", "N", "O", "S")))
  
  if(printplot){
    gp <- ggplot(ft_isodata, aes(x=M_area, y=value)) + 
      geom_point() + 
      geom_smooth(method = "lm") +
      facet_wrap(~name, scales = "free_y") +
      ggtitle(feature_num)
    print(gp)
  }
  
  lmoutput <- split(ft_isodata, ft_isodata$name) %>%
    lapply(lm, formula=value~M_area) %>%
    lapply(summary)
  
  count_ests <- mapply(checkbinom, lminput=lmoutput, n_atoms=c(1,1,1,1,1,2), 
                       prob=c(0.011, 0.00368, 0.00205, 0.0075, 0.0421, 0.011))
  return(count_ests)
}
checkbinom <- function(lminput, n_atoms, prob){
  rsquared <- lminput$r.squared
  if(is.na(rsquared))return(NA)
  if(rsquared<0.99)return(NA)
  coefs <- lminput$coefficients
  norm_factor <- sapply(0:10, dbinom, x=0, prob=prob)
  pred_values <- sapply(0:10, dbinom, x=n_atoms, prob=prob)
  est_slope <- coefs["M_area", "Estimate"]
  (0:10)[which.min(abs(pred_values/norm_factor-est_slope))]
}
rdisop_check <- function(feature_num, final_features, database_formulae){
  feature_msdata <- final_features[final_features$feature==feature_num, ]
  ms1 <- rbind(c(feature_msdata$mzmed, feature_msdata$avgarea),
               c(feature_msdata$mzmed+1.003355, feature_msdata$C13),
               c(feature_msdata$mzmed+1.003355*2, feature_msdata$X2C13),
               c(feature_msdata$mzmed+1.995796, feature_msdata$S34),
               c(feature_msdata$mzmed+0.997035, feature_msdata$N15),
               c(feature_msdata$mzmed+2.004244, feature_msdata$O18))
  ms1 <- ms1[ms1[,2]!=0, , drop=FALSE]
  ms1[,1] <- ms1[,1]-1.007276
  rdoutput <- Rdisop::decomposeIsotopes(masses = ms1[,1], intensities = ms1[,2], 
                                        ppm = ifelse(ms1[1,1]<200, 5*200/ms1[1,1], 5), 
                                        maxisotopes = 2)
  if(is.null(rdoutput)){return(NA)}
  rd_df <- rdoutput %>%
    `[[<-`("isotopes", NULL) %>%
    do.call(what = cbind) %>% 
    as.data.frame(stringsAsFactors=FALSE) %>%
    filter(valid=="Valid"&DBE>=-1) %>%
    filter(formula%chin%database_formulae)
  if(!nrow(rd_df)){return(NA)}
  rd_df %>% slice(1) %>%
    pull(formula) %>%
    unlist()
}




# Grab MSMS data ----
MSMS_files <- "mzMLs/MSMS/" %>%
  list.files(pattern = ".mzML", full.names = TRUE) %>%
  normalizePath() %>%
  `[`(grepl("pos", x = .))
raw_MSMS_data <- lapply(MSMS_files, grabSingleFileMS2) %>%
  mapply(FUN = cbind, as.numeric(gsub(".*DDApos|.mzML", "", MSMS_files)), 
         SIMPLIFY = FALSE) %>%
  lapply(`names<-`, c("rt", "premz", "fragmz", "int", "voltage")) %>%
  do.call(what = rbind) %>% as.data.table()
has_msms <- final_features %>%
  split(.$feature) %>%
  sapply(function(x){
    raw_MSMS_data %>%
      filter(premz%between%pmppm(x$mzmed, ppm = 5)&
               rt%between%(x$rtmed+c(-20, 20))) %>%
      nrow() %>%
      as.logical()
  })
final_features <- mutate(final_features, has_msms=has_msms)
final_features %>% 
  select(c("feature", "mzmed", "rtmed", "avgarea", "has_msms")) %>%
  as.data.frame() %>% 
  head(20)



# Create .mgf files for SIRIUS to read ----
if(dir.exists(sirius_project_dir)){
  unlink(sirius_project_dir, recursive = TRUE)
}
dir.create(sirius_project_dir)
dir.create(paste0(sirius_project_dir, "//raw_files"))
dir.create(paste0(sirius_project_dir, "//output_dir"))

for(feature_num in final_features$feature){
  output_file <- paste0(sirius_project_dir, "\\raw_files\\", feature_num, ".mgf")
  feature_msdata <- final_features[final_features$feature==feature_num, ]
  ms1 <- rbind(c(feature_msdata$mzmed, feature_msdata$avgarea),
               c(feature_msdata$mzmed+1.003355, feature_msdata$C13),
               c(feature_msdata$mzmed+1.003355*2, feature_msdata$X2C13),
               c(feature_msdata$mzmed+1.995796, feature_msdata$S34),
               c(feature_msdata$mzmed+0.997035, feature_msdata$N15),
               c(feature_msdata$mzmed+2.004244, feature_msdata$O18))
  ms1 <- ms1[ms1[,2]!=0, , drop=FALSE]
  ms2 <- raw_MSMS_data[premz%between%pmppm(feature_msdata$mzmed)&
                         rt%between%(feature_msdata$rtmed+c(-10, 10))]
  mgf_maker(feature_msdata = feature_msdata, ms1 = ms1, 
            ms2 = ms2, output_file = output_file)
}



# Run SIRIUS ----
sirius_cmd <- paste0('"C://Program Files//sirius-win64-4.4.17//',
                     'sirius-console-64.exe" ',
                     ' -i "', normalizePath(sirius_project_dir), '//raw_files"',
                     ' -o "', normalizePath(sirius_project_dir), '//output_dir"',
                     ' formula',
                     ' --database PUBCHEM',
                     ' --profile orbitrap',
                     ' --ions-enforced [M+H]+',
                     ' -c 50',
                     ' zodiac',
                     ' fingerid',
                     ' --database bio',
                     ' canopus')
message(sirius_cmd)
system(sirius_cmd)

# Rename "csv"s to actual tsvs to facilitate reading
csv_names <- sirius_project_dir %>%
  list.files(recursive = TRUE) %>%
  grep(pattern = ".csv", value = TRUE) %>%
  paste0(sirius_project_dir, "/", .)
tsv_names <- gsub(pattern = "csv", "tsv", csv_names)
sum(!file.rename(csv_names, tsv_names))

# Read in SIRIUS data ----
sirius_formulas <- read.table(paste0(sirius_project_dir, "/output_dir/formula_",
                                     "identifications.tsv"), sep = "\t",
                              row.names = NULL, header = TRUE, 
                              stringsAsFactors = FALSE) %>%
  `names<-`(c(names(.)[-1], "drop")) %>%
  select(-last_col()) %>%
  mutate(feature=sub(pattern = "_FEATURE1", replacement = "", x = .$id)) %>%
  mutate(feature=sub(pattern = ".*_", replacement = "", x = .$feature)) %>%
  arrange(feature) %>%
  select(feature, sirius_formula=molecularFormula)



# Run Rdisop ----
rdisop_formulas <- final_features$feature %>%
  pbsapply(rdisop_check, final_features = final_features, 
           database_formulae=readRDS("XCMS/data_raw/unique_formulae.rds")) %>%
  unlist() %>%
  cbind(feature=final_features$feature) %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  `names<-`(c("rdisop_formula", "feature"))



# Check isotope matches ----
iso_formulas <- final_features$feature %>%
  pblapply(isocheck, final_peaks=final_peaks) %>%
  do.call(what=rbind) %>% 
  as.data.frame(stringsAsFactors=FALSE) %>% 
  cbind(feature=final_features$feature, stringsAsFactors=FALSE) %>%
  select(paste0(c("C13", "N15", "O18", "S34"), "_area"), "feature") %>%
  `names<-`(c("C", "N", "O", "S", "feature"))

isocheck(feature_num = "FT100", final_peaks = final_peaks, printplot = TRUE)



# Assemble formulas, remove disagreeing ones and duplicates ----
feature_formulas <- final_features %>%
  left_join(sirius_formulas, by="feature") %>%
  left_join(rdisop_formulas, by="feature") %>%
  left_join(iso_formulas, by="feature") %>%
  select(feature, mzmed, rtmed, avgarea, 
         sirius_formula, rdisop_formula, 
         C, N, O, S)

element_formulas <- formula2elements(feature_formulas$sirius_formula)
iso_mismatches <- sapply(c("C", "N", "O", "S"), function(element){
  elem_match <- sapply(seq_along(element_formulas), function(i){
    element_formulas[[i]][element]==feature_formulas[i,][[element]]
  })}) %>% 
  `!`() %>% rowSums(na.rm = TRUE)

best_formulas <- feature_formulas %>%
  filter(rdisop_formula==sirius_formula&!iso_mismatches) %>%
  select(feature, mzmed, rtmed, avgarea, formula=sirius_formula) %>%
  as.data.frame()

duped_features <- lapply(seq_len(nrow(best_formulas)), function(i){
  feature_data <- best_formulas[i,]
  dup_candidates <- best_formulas[
    best_formulas$mzmed%between%pmppm(feature_data$mzmed, ppm = 5)&
      best_formulas$rtmed%between%(feature_data$rtmed+c(-10, 10)),]
  return(dup_candidates[-c(which.max(dup_candidates$avgarea)),])
}) %>% 
  do.call(what = rbind) %>% `[`(duplicated(.),)

best_formulas <- anti_join(x = best_formulas, y=duped_features, by="feature")
write.csv(x = best_formulas, file = "XCMS/data_pretty/feature_formulas.csv")
