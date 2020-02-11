
# Setup things ----
library(pbapply)
library(dplyr)
library(xcms)
setMSnbaseVerbose(TRUE)
register(SerialParam())

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
makeSpectra <- function(eic, filenum){
  scan_list <- split(eic, eic$rt)
  scan_spectra <- lapply(seq_along(scan_list), function(x){
    new("Spectrum1", 
        rt=unique(scan_list[[x]]$rt),
        mz=scan_list[[x]]$mz, 
        intensity=scan_list[[x]]$int, 
        fromFile=as.integer(filenum), 
        acquisitionNum=x,
        scanIndex=x,
        polarity=NA_integer_,
        msLevel=1L)
  })
  spectra_names <- sprintf("%04d", seq_along(scan_list))
  names(scan_spectra) <- paste0("F", sprintf("%02d", unique(eic$file)), ".S", spectra_names)
  scan_spectra
}
makeMSnExp <- function(spectra_list, fnames){
  scans_by_file <- unlist(spectra_list, recursive = FALSE)
  scan_env <- as.environment(unlist(scans_by_file, recursive = FALSE))
  y <- new("MSnExp")
  y@processingData@files <- fnames
  y@assayData <- scan_env
  y@featureData <- AnnotatedDataFrame(
    data.frame(spectrum=1:length(scan_env),
               polarity=rep(1, length(scan_env)),
               totIonCurrent=100000000,
               basePeakIntensity=1000000,
               injectionTime=rnorm(length(scan_env), mean = 50, sd = 5),
               filterString=rep("FTMS + p ESI Full ms [60.0000-900.0000]", length(scan_env)),
               spectrumId=paste0("controllerType=0 controllerNumber=1 scan=", seq_along(scan_env)),
               centroided=rep(TRUE, length(scan_env)),
               scanWindowLowerLimit=rep(100, length(scan_env)),
               scanWindowUpperLimit=rep(101, length(scan_env)),
               row.names = sort(names(scan_env)),
               msLevel = 1L,
               stringsAsFactors = FALSE))
  fData(y)$injectionTime <- 0
  y
}


# Make some data ----
sample_files <- list.files("mzMLs", pattern = "Smp", full.names = TRUE)
ingalls_stans <- read.csv(paste0("https://raw.githubusercontent.com/kheal/Examp",
                                 "le_Untargeted_Metabolomics_Workflow/master/In",
                                 "galls_Lab_Standards.csv"), stringsAsFactors = FALSE) %>%
  filter(Column=="HILIC") %>% filter(ionization_form=="[M+H]") %>% 
  select(Compound.Name, Emperical.Formula, m.z)

raw_data_list <- pblapply(sample_files, grabSingleFileData)
raw_data_list <- lapply(seq_along(raw_data_list), function(x){cbind(file=x, raw_data_list[[x]])})


# Plot it if you like ----
prettyPlotMass <- function(mass=as.numeric(readline(prompt = "Enter a mass:")),
                           eic_list=raw_data_list, linecols=NULL, ppm=4,
                           xlims=NULL, ylims=NULL, maintxt=NULL){
  #Grab the relevant data
  mzr <- pmppm(mass, ppm = ppm)
  mass_eic_list <- lapply(eic_list, function(x){
    filter(x, x$mz>min(mzr)&x$mz<max(mzr))
  })
  if(!any(sapply(mass_eic_list, nrow))){
    plot(1, type="n", main=maintxt)
    text(1, 1, labels = "No data found")
    return()
  }
  mass_eic <- do.call(rbind, mass_eic_list)
  
  #Add some options
  if(!length(linecols)){linecols <- rep("black", length(eic_list))}
  if(!length(xlims)){xlims <- c(min(mass_eic$rt), max(mass_eic$rt))}
  if(!length(ylims)){ylims <- c(0, max(mass_eic$int))}
  if(!length(maintxt)){
    maintxt <- paste(round(min(mass_eic$mz), digits = 5), "-", 
                     round(max(mass_eic$mz), digits = 5))
  }
  
  #Aaaaaand plot
  plot(1, ylim=ylims, xlim=xlims, ylab="Intensity", 
       xlab="Retention time (s)", main=maintxt)
  for(i in seq_along(eic_list)){
    all_rts <- data.frame(rt=unique(eic_list[[i]]$rt))
    mass_subset <- merge.data.frame(mass_eic_list[[i]], all_rts, all.y=TRUE)
    lines(mass_subset$rt, mass_subset$int, type="l",
          ylab="Intensity", xlab="Retention time", col=linecols[i])
  }
}
linecolors <- c("#0000FF33", "#00FF0033")[rep(c(1,1,1,2,2,2), 4)]
prettyPlotMass(mass = as.numeric(ingalls_stans$m.z[ingalls_stans$Compound.Name=="Betaine"]),
               linecols = linecolors)
# for(i in seq_len(nrow(ingalls_stans))){
#   png(filename = paste0("XCMS/standard_chroms/", ingalls_stans$Compound.Name[i], ".png"))
#   prettyPlotMass(mass = as.numeric(ingalls_stans$m.z[i]), linecols = linecolors, 
#                  maintxt=ingalls_stans$Compound.Name[i])
#   dev.off()
# }


# Making my own MSnExp for targeted XCMS things ----
cmpd_master_eic <- do.call(rbind, pblapply(as.numeric(unique(ingalls_stans$m.z)), 
                                           FUN=grabMassData, eic_list=raw_data_list))
file_eics <- split(cmpd_master_eic, cmpd_master_eic$file)
file_spectra <- pblapply(seq_along(file_eics), function(x){
  makeSpectra(file_eics[[x]], filenum = x)
})
new_MSnExp <- makeMSnExp(spectra_list = file_spectra, fnames = basename(sample_files))
writeMSData(new_MSnExp, file = paste0("XCMS/simple_samples/", basename(sample_files)))



# Continuing with targeted analysis ----
sample_files <- list.files("XCMS/simple_samples/", full.names = TRUE)
sample_meta <- gsub(x = basename(sample_files), replacement = "", 
                    pattern = "190715_Smp_FK180310") %>%
  data.frame(sample_name = .,
             sample_group = rep(c(rep("25m", 3), rep("DCM", 3)), 4),
             spindir = rep(c("Cyclonic", "Anti"), each=12),
             stringsAsFactors = FALSE) %>%
  AnnotatedDataFrame()
std_MSnExp <- readMSData(files = list.files("XCMS/simple_samples/", full.names = TRUE), 
                          msLevel. = 1, pdata = sample_meta, mode = "onDisk")
std_MSnExp <- filterRt(std_MSnExp, rt=c(120, 1100))
basic_chrom <- chromatogram(std_MSnExp, aggregationFun = "max")
load("XCMS/base_chrom")
par(mfrow=c(2,1))
par(mar=c(0.1, 2.1, 1.1, 0.1))
plot(base_chrom, col=c("#0000FF33", "#00FF0033")[as.factor(std_MSnExp$sample_group)],
     ylab="", xaxt="n", xlab="", ylim=c(0,100000000))
par(mar=c(2.1, 2.1, 1.1, 0.1))
plot(basic_chrom, col=c("#0000FF33", "#00FF0033")[as.factor(std_MSnExp$sample_group)],
     ylab="", ylim=c(0,100000000))
layout(1)
par(mar=c(2.1, 2.1, 1.1, 1.1))

# Go find that weird peak that's at the surface but not at depth and not in stds ----
library(ggplot2)
all_data <- do.call(rbind, raw_data_list)
filter(all_data, rt>400&rt<450) %>% filter(mz>121.0319&mz<121.0321) %>%
  group_by(file, rt) %>% summarize(bpc=max(int)) %>% 
  mutate(sample_group=c("DCM", "25m")[ceiling(file/3)%%2+1]) %>%
  ggplot(aes(x=rt, y=bpc)) + geom_path(aes(group=file, color=sample_group)) +
  scale_color_manual(values = c("blue", "green"))
filter(all_data, rt>530&rt<560) %>% filter(mz>74.0970&mz<74.0972) %>%
  group_by(file, rt) %>% summarize(bpc=max(int)) %>% 
  mutate(sample_group=c("DCM", "25m")[ceiling(file/3)%%2+1]) %>%
  ggplot(aes(x=rt, y=bpc)) + geom_path(aes(group=file, color=sample_group)) +
  scale_color_manual(values = c("blue", "green"))
filter(all_data, rt>540&rt<570) %>% filter(mz>104.106&mz<104.108) %>%
  group_by(file, rt) %>% summarize(bpc=max(int)) %>% 
  mutate(sample_group=c("DCM", "25m")[ceiling(file/3)%%2+1]) %>%
  ggplot(aes(x=rt, y=bpc)) + geom_path(aes(group=file, color=sample_group)) +
  scale_color_manual(values = c("blue", "green"))
all_data %>% filter(mz>76.075&mz<76.077) %>% filter(rt>275&rt<355) %>%
  group_by(file, rt) %>% summarize(bpc=max(int)) %>% 
  mutate(sample_group=c("DCM", "25m")[ceiling(file/3)%%2+1]) %>%
  ggplot(aes(x=rt, y=bpc)) + geom_path(aes(group=file, color=sample_group)) +
  scale_color_manual(values = c("blue", "green"))
all_data %>% filter(mz>147.075&mz<147.077) %>% filter(rt>580&rt<610) %>%
  group_by(file, rt) %>% summarize(bpc=max(int)) %>% 
  mutate(sample_group=c("DCM", "25m")[ceiling(file/3)%%2+1]) %>%
  ggplot(aes(x=rt, y=bpc)) + geom_path(aes(group=file, color=sample_group)) +
  scale_color_manual(values = c("blue", "green"))

#Generic finder of tall peak in BPC by retention time
rt_of_interest <- 600
filter(all_data, rt<rt_of_interest+5&rt>rt_of_interest-5) %>%
  filter(mz>123.041|mz<123.040) %>%
  group_by(file) %>% summarize(temp=mz[which.max(int)]) %>% 
  group_by(mass=round(temp, digits = 3)) %>% summarize(count=n())

# Spacer ----



