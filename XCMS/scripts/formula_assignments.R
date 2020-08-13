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
final_features %>%
  mutate(has_msms=has_msms) %>%
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
           database_formulae=readRDS("XCMS/unique_formulae.rds")) %>%
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
