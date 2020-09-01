# Script to run SIRIUS (4.4.17), Rdisop, and check isotope signatures and obtain
# high-confidence formula assignments to features
# Expects a feature list (usually created by peakpicking.R)
# Isotope checking expects a peak list (usually created by peakpicking.R)
# Needs bugtesting after WD restructuring


# Grab MSMS data ----
message("Reading in MSMS files...")
raw_MS_data <- pblapply(MSMS_files, grabSingleFileData) %>%
  mapply(FUN = cbind, basename(MSMS_files), SIMPLIFY = FALSE) %>%
  lapply(`names<-`, c("rt", "mz", "int", "file_name")) %>%
  do.call(what = rbind) %>% as.data.table()
raw_MSMS_data <- pblapply(MSMS_files, grabSingleFileMS2) %>%
  mapply(FUN = cbind, basename(MSMS_files), SIMPLIFY = FALSE) %>%
  lapply(`names<-`, c("rt", "premz", "fragmz", "int", "file")) %>%
  do.call(what = rbind) %>% as.data.table() %>%
  mutate(voltage=str_extract(.$file, pattern = "pos\\d+")) %>%
  mutate(voltage=gsub(pattern = "pos", "", .$voltage))


has_msms <- final_features %>%
  split(.$feature) %>%
  pbsapply(function(x){
    raw_MSMS_data[premz%between%pmppm(x$mzmed, ppm = 5)] %>%
      filter(rt%between%(x$rtmed+c(-50, 50))) %>%
      nrow() %>%
      as.logical()
  })



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
sirius_formulas <- sirius_project_dir %>%
  paste0("/output_dir") %>%
  dir(pattern = "formula_candidates", recursive = TRUE, full.names = TRUE) %>%
  pbsapply(read.table, sep = "\t", row.names=NULL, header=TRUE, simplify = FALSE) %>%
  imap(.f = function(x, y){
    cbind(x, feature_num=str_extract(y, "FT\\d+"))
  }) %>%
  lapply(`[[`, "molecularFormula")
names(sirius_formulas) <- str_extract(names(sirius_formulas), "FT\\d+")



# Run Rdisop ----
database_formulae=readRDS("XCMS/unique_formulae.rds")
rdisop_formulas <- final_features$feature %>%
  pbsapply(rdisop_check, final_features = final_features, 
           database_formulae = database_formulae)
rdisop_formulas %>%
  imap(.f = function(x, y){
    cbind(formula=x, feature_num=y)
  }) %>%
  do.call(what=rbind) %>%
  as.data.frame()



# Check isotope matches ----
iso_formulas <- final_features$feature %>%
  pblapply(isocheck, final_peaks=final_peaks) %>%
  do.call(what=rbind) %>% 
  as.data.frame(stringsAsFactors=FALSE) %>% 
  cbind(feature=final_features$feature, stringsAsFactors=FALSE) %>%
  select(paste0(c("C13", "N15", "O18", "S34"), "_area"), "feature") %>%
  `names<-`(c("C", "N", "O", "S", "feature")) %>%
  `rownames<-`(final_features$feature)

isocheck(feature_num = "FT079", final_peaks = final_peaks, printplot = TRUE)



# Assemble formulas, remove disagreeing ones and duplicates ----
inter_formulas <- sapply(names(rdisop_formulas), function(feature_num){
  intersect(rdisop_formulas[[feature_num]], sirius_formulas[[feature_num]])
}, simplify=FALSE)

isochecked_formulas <- lapply(names(inter_formulas), function(feature_num){
  isodata <- iso_formulas[feature_num, ]
  if(!length(inter_formulas[[feature_num]])){
    return(character(0))
  }
  formula_agreements <- inter_formulas[[feature_num]] %>% 
    formula2elements() %>%
    sapply(function(elem_table){
      empty_elements <- c("C", "N", "O", "S")[!c("C", "N", "O", "S")%in%names(elem_table)]
      empty_elements <- `names<-`(numeric(length(empty_elements)), empty_elements)
      elem_table <- c(elem_table, empty_elements)
      c(
        C=elem_table[["C"]]-isodata[["C"]],
        N=elem_table[["N"]]-isodata[["N"]],
        O=elem_table[["O"]]-isodata[["O"]],
        S=elem_table[["S"]]-isodata[["S"]]
      )
  }, simplify = FALSE) %>%
    sapply(sum, na.rm=TRUE) %>%
    sapply(abs)
  best_formula <- inter_formulas[[feature_num]][which(formula_agreements==min(formula_agreements))]
  return(best_formula)
})


best_formulas <- inter_formulas %>%
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

