# Script to de-isotope and de-adduct XCMS-picked peaks
# Called by Control.Rmd

# Find peaks with isotopes ----
# Looks +/- a certain m/z determined by the adduct table
is_peak_iso <- 
  bplapply(split(raw_peaks, raw_peaks$file_name),
           FUN = isIsoAdduct, xdata=xdata_filled,
           grabSingleFileData=grabSingleFileData, checkPeakCor=checkPeakCor, 
           pmppm=pmppm, trapz=trapz, polarity=polarity) %>%
  do.call(what = rbind) %>% as.data.frame()
write.csv(is_peak_iso, file = paste0(intermediate_folder, "is_peak_iso.csv"), row.names = FALSE)
#6.5 minutes

peakshapematch <- is_peak_iso %>%
  group_by(feature) %>%
  summarise(C13_prob=median(C13_match), X2C13_prob=median(X2C13_match), 
            S34_prob=median(S34_match), S33_prob=median(S33_match),
            N15_prob=median(N15_match), O18_prob=median(O18_match),
            Na_prob=median(Na_match), NH4_prob=median(NH4_match), 
            H2O_H_prob=median(H2O_H_match), K_prob=median(K_match),
            X2H_prob=median(X2H_match)) %>%
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

likely_addisos <- peakareamatch$feature[
  which(rowSums(peakareamatch[,names(peakareamatch)!="feature"]>0.95&
                  peakshapematch[,names(peakshapematch)!="feature"]>0.9)>=1)
  ]

addiso_features <- raw_peaks %>%
  group_by(feature) %>%
  summarise(mzmed=median(mz), rtmed=median(rt), avginto=mean(into, na.rm=TRUE)) %>%
  filter(feature%in%likely_addisos)
write.csv(addiso_features, file = paste0(pretty_folder, "addiso_features.csv"), 
          row.names = FALSE)

message(Sys.time()-start_time)
#40 minutes


# Calculate isotopes and adducts for remaining peaks ----
raw_peaks <- read.csv(paste0(intermediate_folder, "raw_peaks.csv"))
addiso_features <- read.csv(paste0(pretty_folder, "addiso_features.csv"))

# For each peak, look for data at +/- each adduct/isotope m/z 
# Also calculate cor while the raw data is being accessed anyway
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
  summarise(C13_cor=median(C13_match), X2C13_cor=median(X2C13_match), 
            S34_cor=median(S34_match), N15_cor=median(N15_match), 
            O18_cor=median(O18_match), S33_cor=median(S33_match),
            K_cor=median(K_match), Na_cor=median(Na_match), 
            NH4_cor=median(NH4_match), H2O_H_cor=median(H2O_H_match), 
            X2H_cor=median(X2H_match)) %>% 
  pivot_longer(cols = -c("feature"), names_to = "addiso", values_to = "cor") %>%
  mutate(addiso=gsub("_cor", "", addiso))
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
  group_by(feature) %>%
  summarize(mzmed=median(mz), rtmed=median(rt), avgarea=mean(M_area)) %>%
  left_join(peak_cors, by="feature") %>%
  left_join(peak_R2s, by=c("feature", "addiso")) %>%
  left_join(peak_slopes, by=c("feature", "addiso")) %>%
  mutate(rel_int=ifelse(cor>0.8&R2>0.9, round(slope*avgarea), 0)) %>%
  select(-c("cor", "R2", "slope")) %>%
  pivot_wider(names_from = addiso, values_from = rel_int)

write.csv(x = complete_peaks, 
          file = paste0(intermediate_folder, "complete_peaks.csv"),
          row.names = FALSE)
write.csv(x = complete_features, 
          file = paste0(intermediate_folder, "complete_features.csv"),
          row.names = FALSE)
message(Sys.time()-start_time)

