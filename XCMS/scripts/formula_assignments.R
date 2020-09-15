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
sirius_cmd <- paste0('"C://Program Files//sirius-win64-4.4.29//',
                     'sirius-console-64.exe" ',
                     ' -i "', normalizePath(sirius_project_dir), '//raw_files"',
                     ' -o "', normalizePath(sirius_project_dir), '//output_dir"',
                     ' formula',
                     ' --database PUBCHEM',
                     ' --profile orbitrap',
                     ' --ions-enforced [M+H]+')
message(sirius_cmd)
system(sirius_cmd)

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
iso_abundance_table <- data.frame(
  isotope=c("C13", "N15", "O18", "X2C13", "S33", "S34"),
  abundance=c(0.011, 0.00368, 0.00205, 0.0075, 0.0421, 0.011),
  n_atoms=c(1,1,1,1,1,2)
)
iso_formulas <- final_features$feature %>% 
  pblapply(isocheck, final_peaks=real_peaks) %>%
  c(list(data.frame(isotope=c("C13", "N15", "O18", "X2C13", "S33", "S34"))), .) %>%
  purrr::reduce(.f=left_join, by="isotope") %>%
  t() %>% as.data.frame() %>% mutate(feature=rownames(.)) %>%
  `colnames<-`(slice(., 1) %>% `[`(-length(.)) %>% c("feature")) %>% 
  slice(-1) %>% select(feature, everything())

isocheck(feature_num = "FT0003", final_peaks = real_peaks, printplot = TRUE)



# Assemble formulas, remove disagreeing ones and duplicates ----
inter_formulas <- sapply(names(rdisop_formulas), function(feature_num){
  intersect(rdisop_formulas[[feature_num]], sirius_formulas[[feature_num]])
}, simplify=FALSE)

isochecked_formulas <- lapply(names(inter_formulas), function(feature_num){
  isodata <- filter(iso_formulas, feature==feature_num)
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
        C=elem_table[["C"]]-as.numeric(isodata[["C13"]]),
        N=elem_table[["N"]]-as.numeric(isodata[["N15"]]),
        O=elem_table[["O"]]-as.numeric(isodata[["O18"]]),
        S=elem_table[["S"]]-as.numeric(isodata[["S34"]])
      )
  }, simplify = FALSE) %>%
    sapply(sum, na.rm=TRUE) %>%
    sapply(abs)
  best_formula <- inter_formulas[[feature_num]][which(formula_agreements==min(formula_agreements))]
  return(best_formula)
}) %>% `names<-`(names(inter_formulas))
  
final_formulas <- data.frame(feature=names(isochecked_formulas), 
           formula=sapply(isochecked_formulas, paste, collapse="; ")) %>%
  left_join(final_features, by="feature") %>%
  select(feature, mzmed, rtmed, avgarea, formula)
