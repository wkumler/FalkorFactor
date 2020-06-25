scrapeMS2 <- function(ms2standata){
  frags <- grep("^[[:digit:]]", ms2standata, value = TRUE)
  frag_consise <- sapply(frags, strsplit, split="\\t", USE.NAMES = FALSE) %>%
    lapply(paste, collapse=":") %>% paste(collapse = "; ")
}
#For positive
pos_stans_ms2_raw <- readLines("Z:/1_QEdata/Will/Ingalls_HILICPos_Standards.txt")
pos_stans_ms2_names <- pos_stans_ms2_raw  %>%
  grep(pattern = "^NAME: ", value = TRUE) %>%
  gsub(pattern = "NAME: Ingalls_|NAME: \\\"Ingalls_|\\\"", replacement = "")
pos_stans_ms2 <- pos_stans_ms2_raw %>%
  split(cumsum(grepl("^NAME:", .))) %>%
  lapply(scrapeMS2) %>%
  unlist() %>% unname() %>%
  cbind(Compound.Name=pos_stans_ms2_names, ms2=., Fraction1="HILICPos") %>%
  as.data.frame() %>%
  filter(!duplicated(Compound.Name))



gp1 <- pos_stans_ms2 %>% 
  filter(grepl(pattern = "^Leucine$", Compound.Name)) %>%
  pull(ms2) %>%
  as.character() %>%
  strsplit(split = "; ") %>%
  data.frame() %>%
  `names<-`("raw") %>%
  separate(col = raw, into = c("mz", "int"), sep = ":") %>%
  mutate(mz=as.numeric(mz)) %>%
  mutate(int=as.numeric(int)) %>%
  ggplot() + geom_segment(aes(x=mz, xend=mz, y=0, yend=int)) +
  xlim(c(50, 150)) + 
  ggtitle("Leucine")

gp2 <- data.frame(mz=c(132.12, 73.08, 60.13, 59.12, 58.11), int=4e+07) %>%
  ggplot() + geom_segment(aes(x=mz, xend=mz, y=0, yend=int)) +
  xlim(c(50, 150)) + 
  ggtitle("TMAP (Pohnert)")

gp3 <- raw_MSMS_data[premz%between%pmppm(131.0946+1.007276, ppm = 5)&
                rt%between%c(280, 330) & voltage==50] %>%
  ggplot() + geom_segment(aes(x=fragmz, xend=fragmz, y=0, yend=int)) +
  xlim(c(50, 150)) + 
  ggtitle("Left peak")

gp4 <- raw_MSMS_data[premz%between%pmppm(131.0946+1.007276, ppm = 5)&
                rt%between%c(470, 500) & voltage==50] %>%
  ggplot() + geom_segment(aes(x=fragmz, xend=fragmz, y=0, yend=int)) +
  xlim(c(50, 150)) + 
  ggtitle("Right peak")

gridExtra::grid.arrange(gp1, gp2, gp3, gp4)


gp <- ggplot() +
  geom_segment(data = data.frame(mz=c(132.12, 73.08, 60.13, 59.12, 58.11), int=4e+07),
               aes(x=mz, xend=mz, y=0, yend=int), color="red") +
  geom_segment(raw_MSMS_data[premz%between%pmppm(131.0946+1.007276, ppm = 5)&
                               rt%between%c(470, 500) & voltage==50],
              mapping = aes(x=fragmz, xend=fragmz, y=0, yend=-int))
ggplotly(gp)

ggplotly(gp1)
ggplotly(gp3)
