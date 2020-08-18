
stan_data <- all_stans[3,]
stan_data <- all_stans %>% filter(mz%between%pmppm(118.0865)) %>% slice(1)

duplicate_masses <- all_stans %>% 
  filter(duplicated(mz)|duplicated(mz, fromLast=TRUE)) %>%
  arrange(as.numeric(mz), as.numeric(rt)) %>%
  select(compound_name, rt, mz)


apply(all_stans, 1, function(stan_data){
  possible_stan_features <- final_features %>%
    filter(mzmed%between%pmppm(as.numeric(stan_data["mz"])))
  stan_vals <- final_peaks %>%
    filter(feature%in%possible_stan_features$feature) %>%
    select(feature, mz, rt, into, sn, file_name, M_area) %>%
    filter(grepl("Std", .$file_name)) %>%
    mutate(stan_mix=str_extract(file_name, "4uMStdsMix[1|2]|H2OinMatrix"))
    
  if(length(unique(stan_vals$stan_mix))<2){
    return(0)
  } else {
    stan_vals %>%
      split(.$feature) %>%
      lapply(aov, formula=M_area~stan_mix) %>%
      lapply(summary) %>%
      lapply(`[[`, 1) %>%
      lapply(`[[`, "Pr(>F)") %>%
      lapply(`[`, 1)
  }
})
