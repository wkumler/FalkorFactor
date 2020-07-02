# Script to find the B-MIS for each standard in the untargeted data set
# Called by Control.Rmd

# all_peaks, falkor_stans, polarity, pmppm all defined in Control.Rmd


# Grab the internal standards and clean up a little
internal_stans <- falkor_stans %>%
  filter(compound_type=="Internal Standard") %>%
  filter(.$polarity==!!polarity) %>% #!! makes sure it's the string being referred to
  mutate(mz=as.numeric(mz)) %>%
  mutate(lower_mz_bound=lapply(mz, pmppm, ppm=10) %>% sapply(`[`, 1)) %>%
  mutate(upper_mz_bound=lapply(mz, pmppm, ppm=10) %>% sapply(`[`, 2)) %>%
  mutate(rt_sec=rt*60)

# For each IS, look in the picked peaks and see if one matches mz & rt
found_stans <- internal_stans %>%
  split(seq_len(nrow(.))) %>%
  lapply(function(i){
    suppressMessages(
      all_peaks %>%
        group_by(feature) %>%
        summarize(mzmed=median(mz), rtmed=median(rt)) %>%
        filter(mzmed%between%c(i$lower_mz_bound,i$upper_mz_bound)) %>%
        mutate(stan=i$compound_name) %>%
        mutate(ppm_diff=(abs(i$mz-.$mzmed)/.$mzmed)*1000000) %>%
        mutate(rt_diff=i$rt_sec-.$rtmed) %>%
        select(feature, stan, mzmed, ppm_diff, rtmed, rt_diff)
    )
  }) %>% do.call(what="rbind") %>% as.data.frame()

# Remove all stans for which more than one peak was found
found_stans <- found_stans[
  !found_stans$stan%in%found_stans$stan[duplicated(found_stans$stan)],]

# Plot it prettily
facet_labels <- found_stans %>%
  split(found_stans$feature) %>%
  sapply(function(i){paste(i$stan, i$feature, sep=": ")})
is_peak_iso %>%
  filter(feature%in%found_stans$feature) %>%
  ggplot() +
  geom_bar(aes(x=file_name, y=M_area), stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
  facet_wrap(~feature, ncol = 1, scales = "free_y",
             labeller = as_labeller(facet_labels))
ggsave(filename = paste0(pretty_folder, "internal_stan_values.pdf"), 
       device = "pdf", height = 15, width = 7.5)



# Step 1: Grab the peak areas from the pooled sample(s) & normalize to injection volume
stan_data <- lapply(found_stans$feature, function(feature_num){
  stan_name <- found_stans[found_stans$feature==feature_num, "stan"]
  is_peak_iso %>%
    filter(feature==feature_num) %>%
    mutate(stan_name=stan_name) %>%
    select(stan_name, feature, mz, rt, M_area, file_name) %>%
    left_join(bionorm_values, by="file_name") %>%
    mutate(bionorm_area=M_area/norm_vol)
})

# Step 2: compare every feature to every standard and calculate min CV
BMIS <- pbsapply(unique(complete_peaks$feature), function(feature_num){
  feature_pooled <- complete_peaks %>%
    filter(feature==feature_num) %>%
    slice(grep(pattern = "Poo", file_name)) %>%
    left_join(bionorm_values, by="file_name") %>%
    mutate(bionorm_area=M_area/norm_vol)
  initial_rsd <- sd(feature_pooled$bionorm_area)/mean(feature_pooled$bionorm_area)
  
  suppressMessages(
    stan_improvements <- stan_data %>%
      do.call(what=rbind) %>%
      slice(grep(pattern = "Poo", file_name)) %>%
      select(file_name, stan_name, stan_bionorm_area=bionorm_area) %>%
      left_join(feature_pooled, by="file_name") %>%
      mutate(stan_norm_area=bionorm_area/stan_bionorm_area) %>%
      group_by(stan_name) %>%
      summarize(MIS_rsd=sd(stan_norm_area, na.rm = TRUE)/
                  mean(stan_norm_area, na.rm=TRUE)) %>%
      ungroup() %>%
      mutate(improvement=(initial_rsd-MIS_rsd)/initial_rsd) %>%
      mutate(initial_rsd=initial_rsd) %>%
      mutate(acceptable=improvement>cut.off) %>%
      rbind(c("None", initial_rsd, 0, initial_rsd, TRUE), .) %>%
      filter(improvement==max(improvement, na.rm = TRUE))
  )
  if(nrow(stan_improvements)>1){
    return("None")
  } else {
    return(stan_improvements$stan_name)
  }
}) %>%
  data.frame(feature=names(.), BMIS=.)

# Step 3: Calculate new peak areas, normalizing to B-MIS
stan_df <- stan_data %>%
  do.call(what = rbind) %>%
  select("BMIS"=stan_name, file_name, bionorm_area) %>%
  rbind(data.frame(BMIS="None", file_name=unique(.$file_name), bionorm_area=1))
BMISed_feature_peaks <- complete_peaks %>%
  left_join(BMIS, by="feature") %>%
  left_join(stan_df, by=c("BMIS", "file_name")) %>%
  arrange(feature) %>%
  group_by(feature) %>%
  mutate(BMISed_area=(M_area/bionorm_area)*mean(bionorm_area, na.rm=TRUE)) %>%
  ungroup() %>%
  select(-bionorm_area) %>%
  filter(!feature%in%found_stans$feature)
BMISed_features <- BMISed_feature_peaks %>%
  group_by(feature) %>%
  summarise(mzmed=median(mz), rtmed=median(rt), BMIS=unique(BMIS), BMIS_avg=mean(BMISed_area, na.rm=TRUE))
write.csv(BMISed_feature_peaks, 
          file = paste0(intermediate_folder, "BMISed_feature_peaks.csv"), 
          row.names = FALSE)



# Write out peak and feature lists ----
write.csv(file = paste0(pretty_folder, "final_features.csv"),
          BMISed_features, row.names = FALSE)
write.csv(file = paste0(pretty_folder, "final_peaks.csv"), 
          BMISed_feature_peaks, row.names = FALSE)
unique(warnings())
message(Sys.time()-start_time)