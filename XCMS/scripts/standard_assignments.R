
stan_guesser <- function(isotope_choice, mix_choice, match_choice, area_choice, rt_choice){
  if(all(is.na(c(isotope_choice, mix_choice, match_choice, area_choice, rt_choice)))){
    return("No peaks found")
  }
  if(!is.na(isotope_choice)){
    return(isotope_choice)
  }
  if(!is.na(match_choice)) {
    return(match_choice)
  }
  if(is.na(mix_choice)){
    mix_options <- NA
  } else if(!nchar(mix_choice)){
    mix_options <- NA
  } else {
    mix_options <- unlist(strsplit(mix_choice, split = "; "))
  }
  if(length(mix_options)==1&!all(is.na(mix_options))){
    return(mix_options)
  }
  
  if(length(mix_options)>1&!all(is.na(mix_options))){
    if(rt_choice==area_choice){
      return(rt_choice)
    }
    if(rt_choice%in%mix_options&!area_choice%in%mix_options){
      return(rt_choice)
    }
    if(area_choice%in%mix_options&!rt_choice%in%mix_options){
      return(rt_choice)
    }
    if(area_choice%in%mix_options&rt_choice%in%mix_options){
      return(paste(c(area_choice, rt_choice), collapse = "; "))
    }
  }
  if(rt_choice==area_choice){
    return(rt_choice)
  }
}

stan_assignments <- all_stans %>%
  mutate(mz=as.numeric(mz)) %>%
  split(.$compound_name) %>%
  pblapply(function(stan_data){
    # dput(stan_data)
    possible_stan_features <- all_features %>%
      filter(mzmed%between%pmppm(as.numeric(stan_data["mz"]))) %>%
      arrange(rtmed)
    if(!nrow(possible_stan_features)){
      return(rep(NA, 5))
    }
    
    # If the compound has an isotopologue, use RT matching
    # If the compound IS an isotopologue, same deal
    isotopologue <- stan_data$compound_name %>%
      paste0(", [15N\\d|13C\\d|2H\\d]") %>%
      grep(x = all_stans$compound_name, value = TRUE)
    if(length(isotopologue)==1){
      isotopo_data <- all_stans %>% 
        filter(compound_name==isotopologue)
      isotopo_features <- all_features %>%
        filter(mzmed%between%pmppm(as.numeric(isotopo_data$mz), 5))
      if(nrow(isotopo_features)==1){
        isotope_choice <- possible_stan_features %>% 
          mutate(rtdiff=abs(rtmed-isotopo_features$rtmed)) %>%
          arrange(rtdiff) %>%
          slice(1) %>%
          pull(feature)
      } else {
        isotope_choice <- NA
      }
    } else if(grepl(pattern = ", 15N|13C|2H", x = stan_data$compound_name)) {
      isotopologue <- gsub(", 15N.*|, 13C.*|, 2H.*", "", stan_data$compound_name)
      isotopo_data <- all_stans %>% 
        filter(compound_name==isotopologue|
                 compound_name==gsub("^D", "", isotopologue)|
                 compound_name==gsub("^L", "", isotopologue))
      isotopo_features <- all_features %>%
        filter(mzmed%between%pmppm(as.numeric(isotopo_data$mz), 5))
      isotope_choice <- possible_stan_features %>% 
        left_join(isotopo_features, by=character()) %>%
        select(-starts_with("mzmed")) %>%
        mutate(rtdiff=abs(rtmed.x-rtmed.y)) %>%
        arrange(rtdiff) %>%
        slice(1) %>%
        pull(feature.x)
    } else {
      isotope_choice <- NA
    }
    
    # If there's a peak that differs between the mixes
    if(is.na(stan_data$mix)){
      mix_choice <- NA
    } else {
      stan_peaks <- all_peaks %>%
        filter(feature%in%possible_stan_features$feature) %>%
        select(feature, mz, rt, into, file_name, M_area) %>%
        filter(grepl("Std", .$file_name)) %>%
        filter(!grepl("H2OinMatrix", .$file_name)) %>%
        mutate(correct_mix=grepl(file_name, pattern = stan_data["mix"])) %>%
        mutate(stan_type=str_extract(file_name, "InH2O|InMatrix"))
      
      # ggplot(stan_peaks) +
      #   geom_boxplot(aes(x=stan_type, y=M_area, color=correct_mix)) +
      #   facet_wrap(~feature, scales = "free_y")
      
      if(!nrow(stan_peaks)){
        return(NA)
      }
      
      mix_peaks <- stan_peaks %>% 
        group_by(feature, correct_mix, stan_type) %>%
        summarise(avgarea=mean(M_area), sdarea=sd(M_area)) %>%
        right_join(expand.grid(
          unique(.$feature),
          unique(.$correct_mix),
          unique(.$stan_type)
        ) %>% `names<-`(c("feature", "correct_mix", "stan_type")),
        by=c("feature", "correct_mix", "stan_type")) %>%
        mutate(avgarea=ifelse(is.na(avgarea), 0, avgarea)) %>%
        mutate(sdarea=ifelse(is.na(sdarea), 0, sdarea))
        
      if(length(unique(mix_peaks$feature))==1){
        mix_choice <- unique(mix_peaks$feature)
      } else {
        mix_choice <- mix_peaks %>%
          split(interaction(.$feature, .$stan_type)) %>%
          lapply(function(v){
            diff <- (v$avgarea[v$correct_mix] - v$avgarea[!v$correct_mix])/mean(v$sdarea)
            data.frame(feature=unique(v$feature), 
                       stan_type=unique(v$stan_type),
                       diff_degree=diff)
          }) %>%
          do.call(what = rbind) %>% `rownames<-`(NULL) %>%
          group_by(feature) %>%
          summarise(correct_mix_peak=mean(diff_degree)) %>%
          filter(correct_mix_peak>10) %>%
          pull(feature) %>%
          paste(collapse = "; ")
      }
    }

    # If there's one peak much closer in RT to expected than the others
    expected_rt <- stan_data$rt*60
    rt_choice <- possible_stan_features %>% 
      mutate(rtdiff=abs(rtmed-expected_rt)) %>%
      arrange(rtdiff) %>%
      slice(1) %>%
      pull(feature)
    
    # If there's same number of features as expected peaks, assume 1:1 and order by RT
    possible_other_stans <- all_stans %>%
      filter(mz%between%pmppm(stan_data$mz, 5)) %>%
      arrange(rt)
    if(nrow(possible_stan_features)==nrow(possible_other_stans)){
      match_choice <- possible_stan_features %>%
        arrange(rtmed) %>%
        cbind(possible_other_stans) %>%
        select(feature, compound_name) %>%
        filter(compound_name==stan_data$compound_name) %>%
        pull(feature)
    } else {
      match_choice <- NA
    }
    
    stan_peaks <- all_peaks %>%
      filter(feature%in%possible_stan_features$feature) %>%
      select(feature, mz, rt, into, file_name, M_area) %>%
      filter(grepl("Std", .$file_name))
    
    if(!nrow(stan_peaks)){
      area_choice <- NA
    } else if(nrow(possible_stan_features)==nrow(possible_other_stans)&
              nrow(possible_stan_features)>1) {
      area_choice <- NA
    } else {
      area_choice <- stan_peaks %>% 
        group_by(feature) %>%
        summarise(avgarea=mean(M_area)) %>%
        arrange(desc(avgarea)) %>%
        slice(1) %>%
        pull(feature)
    }
    
    best_guess <- stan_guesser(isotope_choice, mix_choice, match_choice, area_choice, rt_choice)
    
    return(data.frame(
      best_guess=best_guess,
      isotope_validated=isotope_choice,
      rt_matchup=match_choice,
      mix_matched=mix_choice,
      closer_rt=rt_choice,
      area_choice=area_choice
    ))
  }) %>%
  do.call(what=rbind) %>%
  mutate(compound_name=rownames(.)) %>%
  select(compound_name, everything())


stan_assignments %>%
  filter(!is.na(.$isotope_validated))

stan_assignments[(duplicated(stan_assignments$best_guess, fromLast=TRUE)|
                   duplicated(stan_assignments$best_guess))&
                   !is.na(stan_assignments$best_guess), ] %>%
  split(.$best_guess)

stan_assignments %>%
  mutate(compound_name=rownames(.)) %>%
  `rownames<-`(NULL) %>%
  left_join(all_stans) %>%
  arrange(as.numeric(mz), as.numeric(rt)) %>%
  select(compound_name, best_guess, 1:6)
