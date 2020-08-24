
duplicate_masses <- all_stans %>% 
  filter(duplicated(mz)|duplicated(mz, fromLast=TRUE)) %>%
  arrange(as.numeric(mz), as.numeric(rt)) %>%
  select(compound_name, rt, mz)


stan_assignments <- pbapply(all_stans, 1, function(stan_data){
  possible_stan_features <- final_features %>%
    filter(mzmed%between%pmppm(as.numeric(stan_data["mz"])))
  stan_vals <- final_peaks %>%
    filter(feature%in%possible_stan_features$feature) %>%
    select(feature, mz, rt, into, sn, file_name, M_area) %>%
    filter(grepl("Std", .$file_name)) %>%
    mutate(correct_mix=grepl(file_name, pattern = stan_data["mix"])) %>%
    split(.$feature)
  
  stan_diffs <- sapply(stan_vals, function(x){
    if(sum(x$correct_mix)<2|sum(!x$correct_mix)<2){
      return(0)
    } else {
      return(t.test(M_area~correct_mix, data = x)$p.value)
    }
  })
})

all_stans[which(sapply(stan_assignments, length)==5),]



# Glycine (5 possible peaks)
stan_row <- 47
stan_data <- all_stans[stan_row,]
filter(all_stans, mz%between%pmppm(as.numeric(stan_data["mz"])))
possible_stan_features <- final_features %>%
  filter(mzmed%between%pmppm(as.numeric(stan_data["mz"]))) %>%
  arrange(rtmed)
stan_vals <- final_peaks %>%
  filter(feature%in%possible_stan_features$feature) %>%
  select(feature, mz, rt, into, sn, file_name, M_area) %>%
  filter(grepl("Std", .$file_name)) %>%
  mutate(correct_mix=grepl(file_name, pattern = stan_data["mix"])) %>%
  split(.$feature)
stan_assignments[[stan_row]] %>% data.frame(feature=names(.), pval=.) %>% 
  left_join(possible_stan_features) %>% arrange(rtmed)
read.csv("C:/Users/willi/Downloads/76.03994.csv") %>%
  mutate(fillcol=str_extract(fileid, "Blk|Smp|Std")) %>%
  filter(fillcol=="Std") %>%
  mutate(stdtype=str_extract(fileid, "Mix[1|2]|H2O")) %>%
  ggplot() + 
  geom_polygon(aes(x=rt, y=int, fill=stdtype)) +
  geom_line(aes(x=rt, y=int, group=stdtype)) +
  scale_y_continuous(oob = scales::rescale_none, limits = c(0, 500000)) +
  xlim(c(100, 750)) + 
  facet_wrap(~fileid) +
  scale_fill_viridis_d(alpha = 0.7, option = "C")




all_stans[which(sapply(stan_assignments, length)==4),]

# Allopurinol, hypoxanthine (4 possible peaks)
all_stans[c(91, 106),]
stan_row <- 91
stan_data <- all_stans[stan_row,]
possible_stan_features <- final_features %>%
  filter(mzmed%between%pmppm(as.numeric(stan_data["mz"])))
stan_vals <- final_peaks %>%
  filter(feature%in%possible_stan_features$feature) %>%
  select(feature, mz, rt, into, sn, file_name, M_area) %>%
  filter(grepl("Std", .$file_name)) %>%
  mutate(correct_mix=grepl(file_name, pattern = stan_data["mix"])) %>%
  split(.$feature)
stan_assignments[[stan_row]] %>% data.frame(feature=names(.), pval=.) %>% 
  left_join(possible_stan_features) %>% arrange(rtmed)
read.csv("C:/Users/willi/Downloads/137.046336.csv") %>%
  mutate(fillcol=str_extract(fileid, "Blk|Smp|Std")) %>%
  filter(fillcol=="Std") %>%
  mutate(stdtype=str_extract(fileid, "Mix[1|2]|H2O")) %>%
  ggplot() + 
  geom_polygon(aes(x=rt, y=int, fill=stdtype)) +
  geom_line(aes(x=rt, y=int, group=stdtype)) +
  facet_wrap(~fileid) +
  scale_y_continuous(oob = scales::rescale_none, limits = c(0, 2000000)) +
  xlim(200, 500) +
  scale_fill_viridis_d(alpha = 0.7, option = "C")

# Betonicine/turicine
all_stans[c(9, 85),]
stan_row <- 9
stan_data <- all_stans[stan_row,]
filter(all_stans, mz%between%pmppm(as.numeric(stan_data["mz"])))
possible_stan_features <- final_features %>%
  filter(mzmed%between%pmppm(as.numeric(stan_data["mz"]))) %>%
  arrange(rtmed)
stan_vals <- final_peaks %>%
  filter(feature%in%possible_stan_features$feature) %>%
  select(feature, mz, rt, into, sn, file_name, M_area) %>%
  filter(grepl("Std", .$file_name)) %>%
  mutate(correct_mix=grepl(file_name, pattern = stan_data["mix"])) %>%
  split(.$feature)
stan_assignments[[stan_row]] %>% data.frame(feature=names(.), pval=.) %>% 
  left_join(possible_stan_features) %>% arrange(rtmed)
read.csv("C:/Users/willi/Downloads/160.0969.csv") %>%
  mutate(fillcol=str_extract(fileid, "Blk|Smp|Std")) %>%
  filter(fillcol=="Std") %>%
  mutate(stdtype=str_extract(fileid, "Mix[1|2]|H2O")) %>%
  ggplot() + 
  geom_polygon(aes(x=rt, y=int, fill=stdtype)) +
  geom_line(aes(x=rt, y=int, group=stdtype)) +
  scale_y_continuous(oob = scales::rescale_none, limits = c(0, 1000000)) +
  xlim(300, 600) +
  facet_wrap(~fileid) +
  scale_fill_viridis_d(alpha = 0.7, option = "C")


# Leucine, isoleucine, TMAP
all_stans[c(55, 56, 98),]
stan_row <- 55
stan_data <- all_stans[stan_row,]
filter(all_stans, mz%between%pmppm(as.numeric(stan_data["mz"]))) %>%
  arrange(rt)
possible_stan_features <- final_features %>%
  filter(mzmed%between%pmppm(as.numeric(stan_data["mz"]))) %>%
  arrange(rtmed)
stan_vals <- final_peaks %>%
  filter(feature%in%possible_stan_features$feature) %>%
  select(feature, mz, rt, into, sn, file_name, M_area) %>%
  filter(grepl("Std", .$file_name)) %>%
  mutate(correct_mix=grepl(file_name, pattern = stan_data["mix"])) %>%
  split(.$feature)
stan_assignments[[stan_row]] %>% data.frame(feature=names(.), pval=.) %>% 
  left_join(possible_stan_features) %>% arrange(rtmed)
read.csv("C:/Users/willi/Downloads/132.102.csv") %>%
  mutate(fillcol=str_extract(fileid, "Blk|Smp|Std")) %>%
  filter(fillcol=="Std") %>%
  mutate(stdtype=str_extract(fileid, "Mix[1|2]|H2O")) %>%
  ggplot() + 
  geom_polygon(aes(x=rt, y=int, fill=stdtype)) +
  geom_line(aes(x=rt, y=int, group=stdtype)) +
  xlim(300, 600) + 
  scale_y_continuous(oob = scales::rescale_none, limits = c(0, 10000000)) +
  facet_wrap(~fileid) +
  scale_fill_viridis_d(alpha = 0.7, option = "C")




all_stans[which(sapply(stan_assignments, length)==3),]

# Sarcosine, alanine, beta-alanine
all_stans[c(18, 69, 93),]
stan_row <- 18
stan_data <- all_stans[stan_row,]
filter(all_stans, mz%between%pmppm(as.numeric(stan_data["mz"]))) %>%
  arrange(rt)
possible_stan_features <- final_features %>%
  filter(mzmed%between%pmppm(as.numeric(stan_data["mz"]))) %>%
  arrange(rtmed)
stan_vals <- final_peaks %>%
  filter(feature%in%possible_stan_features$feature) %>%
  select(feature, mz, rt, into, sn, file_name, M_area) %>%
  filter(grepl("Std", .$file_name)) %>%
  mutate(correct_mix=grepl(file_name, pattern = stan_data["mix"])) %>%
  split(.$feature)
stan_assignments[[stan_row]] %>% data.frame(feature=names(.), pval=.) %>% 
  left_join(possible_stan_features) %>% arrange(rtmed)
read.csv("C:/Users/willi/Downloads/90.05551.csv") %>%
  mutate(fillcol=str_extract(fileid, "Blk|Smp|Std")) %>%
  filter(fillcol=="Std") %>%
  mutate(stdtype=str_extract(fileid, "Mix[1|2]|H2O")) %>%
  ggplot() + 
  geom_polygon(aes(x=rt, y=int, fill=stdtype)) +
  geom_line(aes(x=rt, y=int, group=stdtype)) +
  xlim(500, 700) +
  scale_y_continuous(oob = scales::rescale_none, limits = c(0, 10000000)) +
  facet_wrap(~fileid) +
  scale_fill_viridis_d(alpha = 0.7, option = "C")

# Cytosine
all_stans[36,]
stan_row <- 36
stan_data <- all_stans[stan_row,]
filter(all_stans, mz%between%pmppm(as.numeric(stan_data["mz"]))) %>%
  arrange(rt)
possible_stan_features <- final_features %>%
  filter(mzmed%between%pmppm(as.numeric(stan_data["mz"]))) %>%
  arrange(rtmed)
stan_vals <- final_peaks %>%
  filter(feature%in%possible_stan_features$feature) %>%
  select(feature, mz, rt, into, sn, file_name, M_area) %>%
  filter(grepl("Std", .$file_name)) %>%
  mutate(correct_mix=grepl(file_name, pattern = stan_data["mix"])) %>%
  split(.$feature)
stan_assignments[[stan_row]] %>% data.frame(feature=names(.), pval=.) %>% 
  left_join(possible_stan_features) %>% arrange(rtmed)
read.csv("C:/Users/willi/Downloads/112.051087.csv") %>%
  mutate(fillcol=str_extract(fileid, "Blk|Smp|Std")) %>%
  filter(fillcol=="Std") %>%
  mutate(stdtype=str_extract(fileid, "Mix[1|2]|H2O")) %>%
  ggplot() + 
  geom_polygon(aes(x=rt, y=int, fill=stdtype)) +
  geom_line(aes(x=rt, y=int, group=stdtype)) +
  xlim(300, 500) +
  scale_y_continuous(oob = scales::rescale_none, limits = c(0, 10000000)) +
  facet_wrap(~fileid) +
  scale_fill_viridis_d(alpha = 0.7, option = "C")

# Proline
all_stans[66,]
stan_row <- 66
stan_data <- all_stans[stan_row,]
filter(all_stans, mz%between%pmppm(as.numeric(stan_data["mz"]))) %>%
  arrange(rt)
possible_stan_features <- final_features %>%
  filter(mzmed%between%pmppm(as.numeric(stan_data["mz"]))) %>%
  arrange(rtmed)
stan_vals <- final_peaks %>%
  filter(feature%in%possible_stan_features$feature) %>%
  select(feature, mz, rt, into, sn, file_name, M_area) %>%
  filter(grepl("Std", .$file_name)) %>%
  mutate(correct_mix=grepl(file_name, pattern = stan_data["mix"])) %>%
  split(.$feature)
stan_assignments[[stan_row]] %>% data.frame(feature=names(.), pval=.) %>% 
  left_join(possible_stan_features) %>% arrange(rtmed)
read.csv("C:/Users/willi/Downloads/116.071154.csv") %>%
  mutate(fillcol=str_extract(fileid, "Blk|Smp|Std")) %>%
  filter(fillcol=="Std") %>%
  mutate(stdtype=str_extract(fileid, "Mix[1|2]|H2O")) %>%
  ggplot() + 
  geom_polygon(aes(x=rt, y=int, fill=stdtype)) +
  geom_line(aes(x=rt, y=int, group=stdtype)) +
  scale_y_continuous(oob = scales::rescale_none, limits = c(0, 1000000)) +
  facet_wrap(~fileid) +
  scale_fill_viridis_d(alpha = 0.7, option = "C")
