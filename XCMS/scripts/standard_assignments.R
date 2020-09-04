
duplicate_masses <- all_stans %>% 
  filter(duplicated(mz)|duplicated(mz, fromLast=TRUE)) %>%
  arrange(as.numeric(mz), as.numeric(rt)) %>%
  select(compound_name, rt, mz)

stan_assignments <- all_stans %>%
  split(.$compound_name) %>%
  pblapply(function(stan_data){
    dput(stan_data)
    possible_stan_features <- all_features %>%
      filter(mzmed%between%pmppm(as.numeric(stan_data["mz"]))) %>%
      arrange(rtmed)
    
    # If the compound has an isotopologue, use RT matching
    isotopologue <- stan_data$compound_name %>%
      paste0(", [15N|13C|2H]") %>%
      grep(x = all_stans$compound_name, value = TRUE)
    if(length(isotopologue)==1){
      isotopo_data <- all_stans %>% 
        filter(compound_name==isotopologue)
      isotopo_features <- all_features %>%
        filter(mzmed%between%pmppm(isotopo_data$mz, 5))
      if(nrow(isotopo_features)==1){
        isotope_choice <- possible_stan_features %>% 
          mutate(rtdiff=abs(rtmed-isotopo_features$rtmed)) %>%
          arrange(desc(rtdiff)) %>%
          slice(1) %>%
          pull(feature)
      } else {
        isotope_choice <- NA
      }
    } else {
      isotope_choice <- NA
    }
    
    # If there's a peak that differs between the mixes
    if(is.na(stan_data$mix)){
      mix_choice <- NA
    } else {
      stan_pvals <- real_peaks %>%
        filter(feature%in%possible_stan_features$feature) %>%
        select(feature, mz, rt, into, file_name, M_area) %>%
        filter(grepl("Std", .$file_name)) %>%
        mutate(correct_mix=grepl(file_name, pattern = stan_data["mix"])) %>%
        mutate(stan_type=str_extract(file_name, "InH2O|InMatrix|H2OinMatrix")) %>%
        split(.$feature) %>%
        sapply(function(x){
          if(sum(x$correct_mix)<2|sum(!x$correct_mix)<2){
            return(0)
          } else {
            aov_res <- summary(aov(x, formula = M_area~correct_mix+stan_type))
            pvals <- aov_res[[1]]$`Pr(>F)`
            names(pvals) <- c("correct_mix", "stan_type", "Residuals")
            pvals
          }
        }, simplify = FALSE)
      mix_choice <- names(stan_pvals)[
        which(sapply(stan_pvals, function(x)x["correct_mix"]<0.05&x["stan_type"]>0.05))
      ]
    }

    # If there's one peak much closer in RT to expected than the others
    expected_rt <- stan_data$rt*60
    rt_choice <- possible_stan_features %>% 
      mutate(rtdiff=abs(rtmed-expected_rt)) %>%
      arrange(rtdiff) %>%
      slice(1) %>%
      pull(feature)
  })
