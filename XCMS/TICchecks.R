# How much variation in the TIC have we captured in our standards?
library(ggplot2)

# Functions ----
pmppm <- function(mass, ppm=4){c(mass*(1-ppm/1000000), mass*(1+ppm/1000000))}
grabSingleFileData <- function(filename){
  msdata <- mzR:::openMSfile(filename)
  fullhd <- mzR::header(msdata)
  spectra_list <- lapply(seq_len(nrow(fullhd)), function(x){
    given_peaks <- mzR::peaks(msdata, x)
    rtime <- fullhd[x, "retentionTime"]
    return(cbind(rtime, given_peaks))
  })
  all_data <- `names<-`(as.data.frame(do.call(rbind, spectra_list)), 
                        c("rt", "mz", "int"))
  return(all_data)
}
grabMassData <- function(mass, eic_list=raw_data_list, ppm=4){
  mzr <- pmppm(mass, ppm = ppm)
  mass_eic_list <- lapply(eic_list, function(x){
    filter(x, x$mz>min(mzr)&x$mz<max(mzr))
  })
  do.call(rbind, mass_eic_list)
}



# Read in the raw data ----
sample_files <- list.files("mzMLs", pattern = "Smp", full.names = TRUE)
raw_data_list <- pblapply(sample_files, grabSingleFileData)
raw_data_list <- lapply(seq_along(raw_data_list), function(x){cbind(file=x, raw_data_list[[x]])})
all_data <- do.call(rbind, raw_data_list) %>% filter(rt>120&rt<1100)

# Grab the standards & create stan_data
ingalls_stans <- read.csv(paste0("https://raw.githubusercontent.com/kheal/Examp",
                                 "le_Untargeted_Metabolomics_Workflow/master/In",
                                 "galls_Lab_Standards.csv"), stringsAsFactors = FALSE) %>%
  filter(Column=="HILIC") %>% filter(ionization_form=="[M+H]") %>% 
  select(Compound.Name, Emperical.Formula, m.z)
stan_masses <- as.numeric(ingalls_stans$m.z)

shitpeaks <- all_data %>% filter(rt>900&rt<901) %>% 
  mutate(mz=round(mz, digits = 3)) %>%
  group_by(mz) %>% summarize(TIC=sum(int)) %>% 
  arrange(desc(TIC)) %>% slice(1:100) %>% pull(mz)
for(i in shitpeaks$mz){
  gp <- all_data %>% filter(mz>min(pmppm(i))&mz<max(pmppm(i))) %>% 
    group_by(file, rt) %>% summarize(TIC=sum(int)) %>%
    mutate(sample_group=c("DCM", "25m")[ceiling(file/3)%%2+1]) %>%
    ggplot(aes(x=rt, y=TIC, group=file, color=sample_group)) + 
    geom_line() + scale_color_manual(values = c("#0000FF55", "#00FF0055")) +
    theme_bw() + ggtitle(i)
  print(gp)
  readline(prompt = "Press")
}
to_eliminate <- unique(c(stan_masses, shitpeaks))

cleaned_data <- all_data
for(i in seq_along(to_eliminate)){
  massmaxmin <- pmppm(to_eliminate[i])
  cleaned_data <- filter(cleaned_data, mz<min(massmaxmin)|mz>max(massmaxmin))
  gp <- cleaned_data %>% group_by(file, rt) %>% summarize(TIC=sum(int)) %>%
    mutate(sample_group=c("DCM", "25m")[ceiling(file/3)%%2+1]) %>%
    ggplot(aes(x=rt, y=TIC, group=file, color=sample_group)) + 
    geom_line() + scale_color_manual(values = c("#0000FF55", "#00FF0055")) +
    theme_bw()
  ggsave(filename = paste0("XCMS/removalfigs/", sprintf("%03d", i), "mass_", 
                           gsub("\\.", "_", round(to_eliminate[i], digits=3)), ".png"))
}


# Plot the data ----
all_data %>% group_by(file, rt) %>% summarize(TIC=sum(int)) %>%
  mutate(sample_group=c("DCM", "25m")[ceiling(file/3)%%2+1]) %>%
  ggplot(aes(x=rt, y=TIC, group=file, color=sample_group)) + 
  geom_line() + scale_color_manual(values = c("#0000FF55", "#00FF0055")) +
  theme_bw()
stan_data %>% group_by(file, rt) %>% summarize(TIC=sum(int)) %>%
  mutate(sample_group=c("DCM", "25m")[ceiling(file/3)%%2+1]) %>%
  ggplot(aes(x=rt, y=TIC, group=file, color=sample_group)) + 
  geom_line() + scale_color_manual(values = c("#0000FF55", "#00FF0055")) +
  theme_bw()
cleaned_data %>% group_by(file, rt) %>% summarize(TIC=sum(int)) %>%
  mutate(sample_group=c("DCM", "25m")[ceiling(file/3)%%2+1]) %>%
  ggplot(aes(x=rt, y=TIC, group=file, color=sample_group)) + 
  geom_line() + scale_color_manual(values = c("#0000FF55", "#00FF0055")) +
  theme_bw()
