

apply(all_stans, 1, function(stan_data){
  possible_stan_features <- final_features %>%
    filter(mzmed%between%pmppm(as.numeric(stan_data["mz"])))
  final_peaks %>%
    filter(feature%in%possible_stan_features$feature) %>%
    select(feature, mz, rt, into, sn, file_name, M_area) %>%
    group_by(feature) %>%
    filter(grepl("Std", .$file_name)) %>%
    mutate(stan_mix=str_extract(file_name, "4uMStdsMix[1|2]|H2OinMatrix")) %>%
    summarize(p_anova=summary(aov(formula = M_area~stan_mix, data = .))[[1]]$`Pr(>F)`[1])
})
