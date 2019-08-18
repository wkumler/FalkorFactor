# Script to process Ingalls Lab Standards to a friendlier format

# Setup things ----

library(dplyr)
library(Rdisop)

# Grab the original file - not sure where from originally, but currently safe
# in UW folder.
stds <- read.csv("Ingalls_Lab_Standards.csv", stringsAsFactors = F)

getM <- function(mol){
  return(getMolecule(mol)$exactmass)
}

# Filter out unhelpful information and convert rt to seconds

HILIC_pos_stds_mix1 <- stds %>%
  filter(Fraction1=="HILICPos") %>%
  filter(HILICMix=="Mix1") %>%
  filter(ionization_form!="[M]") %>%
  mutate("rt.sec"=RT..min.*60) %>%
  select(Compound.Type, Compound.Name, Emperical.Formula, RT..min., rt.sec, m.z) %>%
  mutate(mass=sapply(as.character(Emperical.Formula), getM)) %>%
  mutate(premass=m.z-mass) %>%
  arrange(mass)

write.csv(HILIC_pos_stds_mix1, file = "Ingalls_Lab_Standards_Will.csv")


# Get rid of compounds not found in the raw data 
# Needs to be checked on MSDIAL raw data

# clean_pos_stds <- HILIC_pos_stds_mix1 %>%
#   filter(!Compound.Name%in%c("Pyridinedicarboxylic acid", "Succinylglycine", 
#                              "Acetylglutamic acid", "Ergothioneine",
#                              "Glutamylphenylalanine"))
# write.csv(clean_pos_stds, file = "Ingalls_Lab_Standards_Will_clean.csv")
