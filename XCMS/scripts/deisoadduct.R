# Script to de-isotope and de-adduct XCMS-picked peaks
# Called by Control.Rmd

# raw_peaks, xdata_filled, area_remove_threshold, shape_remove_threshold
# area_find_threshold, shape_find_threshold all come from Control.Rmd

# Find peaks likely to be isotopes ----
# Looks +/- a certain m/z determined by the adduct table

message("Identifying isotope and adduct peaks... ")
start_time <- Sys.time()
is_peak_iso <- 
  bplapply(split(raw_peaks, raw_peaks$file_name),
           FUN = isIsoAdduct, xdata=xdata_filled,
           grabSingleFileData=grabSingleFileData, checkPeakCor=checkPeakCor, 
           pmppm=pmppm, trapz=trapz, polarity=polarity) %>%
  do.call(what = rbind) %>% as.data.frame()

peakshapematch <- is_peak_iso %>%
  group_by(feature) %>%
  summarise(across(contains("match"), median)) %>%
  as.data.frame()

peakareamatch <- lapply(unique(is_peak_iso$feature), function(i){
  feature_areas <- is_peak_iso[is_peak_iso$feature==i,]
  area_cols <- grep(pattern = "area$", names(feature_areas), value = TRUE)[-1]
  sapply(area_cols, function(x){
    suppressWarnings(cor(feature_areas$M_area, feature_areas[[x]]))
  })
}) %>% 
  do.call(what=rbind) %>% 
  `[<-`(is.na(.), 0) %>% 
  as.data.frame(stringsAsFactors=FALSE) %>%
  mutate(feature=unique(is_peak_iso$feature)) %>%
  select(feature, everything()) %>%
  arrange(feature)

likely_addisos <- peakareamatch[
  which(rowSums(peakareamatch[,names(peakareamatch)!="feature"]>
                  area_remove_threshold&
                  peakshapematch[,names(peakshapematch)!="feature"]>
                  shape_remove_threshold)>=1),
  ]
likely_addisos$adduct_type <- likely_addisos %>%
  split(seq_len(nrow(.))) %>%
  sapply(function(i){
  names(which.max(i[-1]))
}, USE.NAMES = FALSE) %>%
  gsub(pattern = "_area", replacement = "")

addiso_peaks <- raw_peaks %>%
  left_join(likely_addisos[,c("feature", "adduct_type")], by = "feature") %>%
  left_join(is_peak_iso %>% select(feature, file_name, M_area), 
            by = c("feature", "file_name")) %>%
  filter(!is.na(adduct_type)) %>%
  as.data.frame()

message("Time to remove isotopes and adducts: ", 
        round(Sys.time()-start_time, digits = 2), " min")



# Find the adduct/isotope data for all real peaks ----
# For each peak, look for data at +/- each adduct/isotope m/z 
# Also calculate cor while the raw data is being accessed anyway

message("Finding isotopes and adducts for remaining peaks... ")
start_time <- Sys.time()
complete_peaks <- raw_peaks %>%
  filter(!feature%in%addiso_features$feature) %>%
  split(.$file_name) %>%
  bplapply(FUN = findIsoAdduct, xdata=xdata_filled,
           grabSingleFileData=grabSingleFileData, checkPeakCor=checkPeakCor, 
           pmppm=pmppm, trapz=trapz, polarity=polarity) %>%
  do.call(what = rbind) %>% as.data.frame() %>% 
  `rownames<-`(NULL) %>% arrange(feature)

# Calculate median cor for each FEATURE from the various peak cors
peak_cors <- complete_peaks %>%
  group_by(feature) %>%
  summarise(across(contains("match"), median)) %>%
  pivot_longer(cols = -c("feature"), names_to = "addiso", values_to = "cor") %>%
  mutate(addiso=gsub("_match", "", addiso))
# For each feature, plot adduct/iso areas against OG peak areas
# Run lm() to get best fit line slope and R-squared
peak_slope_R2 <- lapply(unique(complete_peaks$feature), function(i){
  feature_areas <- complete_peaks[complete_peaks$feature==i,]
  area_cols <- grep(pattern = "area$", names(feature_areas), value = TRUE)[-1]
  area_outputs <- lapply(area_cols, function(x){
    lmoutput <- lm(feature_areas[[x]]~feature_areas$M_area)
    useful_info <- c(r2=summary(lmoutput)$r.squared, 
                     slope=lmoutput$coefficients["feature_areas$M_area"])
    return(useful_info)
  })
  lapply(area_outputs, c) %>%
    unlist() %>%
    `names<-`(paste0(rep(gsub("area", "", area_cols), each=2), c("R2", "slope")))
}) %>% 
  do.call(what=rbind) %>% `[<-`(is.na(.), 0) %>% 
  as.data.frame(stringsAsFactors=FALSE) %>%
  mutate(feature=unique(complete_peaks$feature)) %>%
  select(feature, everything()) %>%
  arrange(feature)
# Separate out R-squareds and slopes (easier to do here than after merging)
peak_R2s <- peak_slope_R2 %>%
  select(1, grep("R2", names(peak_slope_R2))) %>%
  pivot_longer(cols = -c("feature"), names_to = "addiso", values_to = "R2") %>%
  mutate(addiso=gsub("_R2", "", addiso))
peak_slopes <- peak_slope_R2 %>%
  select(1, grep("slope", names(peak_slope_R2))) %>%
  pivot_longer(cols = -c("feature"), names_to = "addiso", values_to = "slope") %>%
  mutate(addiso=gsub("_slope", "", addiso))


# Establish thresholds for "yes, this is probably an adduct"
# If above threshold, return peak area as relative intensity
# If below, return nothing
# Essentially produces a cleaned up MS1 spectrum with only adducts/isotopes
complete_features <- complete_peaks %>% 
  left_join(peak_cors, by="feature") %>%
  left_join(peak_R2s, by=c("feature", "addiso")) %>%
  left_join(peak_slopes, by=c("feature", "addiso")) %>%
  mutate(rel_int=ifelse(cor>shape_find_threshold&
                          R2>area_find_threshold, 
                        round(slope*M_area), 0)) %>%
  select(-c("cor", "R2", "slope")) %>%
  pivot_wider(names_from = addiso, values_from = rel_int)

message("Time to find new isotopes and adducts: ", 
        round(Sys.time()-start_time, digits = 2), " min")
