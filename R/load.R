
pacman::p_load(car, dplyr, tidyr, nlme, nlshelper, tibble, magicaxis, plantecophys)



if(!dir.exists("download"))dir.create("download")
if(!dir.exists("output"))dir.create("output")
source("R/functions.R")
source("R/figures.R")

# Lin et al. 2015
linfn <- "download/lin2015.csv"
if(!file.exists(linfn)){
  download.file("https://ndownloader.figshare.com/files/1886204", linfn)
}

lin2015 <- read.csv(linfn, stringsAsFactors = FALSE) %>%
  mutate(fitgroup = paste(Datacontrib, Species, sep="_"),
         Ci = CO2S - Photo/(Cond/1.6),
         CiCa = Ci / CO2S,
         BBopti = Photo/(CO2S*sqrt(VPD)))

lin2015a <- group_by(lin2015, fitgroup) %>%
  filter(n() > 15,
         PARin < 2500,
         PARin >= 0,
         Datacontrib != "Derek Eamus",  # sorry derek you are an outlier
         fitgroup != "Jon Bennie_Zornia")

# Lin 2015 fits
lin2015coef <- fits_lin2015(lin2015a)


# Kerstiens 1996
kerst <- read.csv("data/kerstiens1996_table1.csv", stringsAsFactors = FALSE) %>%
  filter(gmin > 0,
         method %in% c("A","C","D")) %>%
  mutate(method = dplyr::recode(method, A="gcut_isol", C="gcut_seal", D="gmin"),
         gmin = 10^5 * 10^-3 * gmin / 40)   # convert to mmol m-2 s-1 assume 40 mol m-3

# Miner et al. 2016
miner <- read.csv("data/Miner_table1.csv")


# Lombardozzi et al 2017
lombar <- read.csv("data/lombardozzi_gnight.csv") %>%
  rename(gmin = gnight) %>%
  filter(gmin > 0) %>%
  mutate(method="gnight")


g0s <- data.frame(gmin=1000 * c(lin2015coef$g0, miner$g0),
                  method="g0", stringsAsFactors = FALSE) %>%
  filter(!is.na(gmin),
         gmin > 0)

gdfr <- bind_rows(kerst, lombar, g0s) %>%
  mutate(method = factor(method, levels=c("gcut_isol","gcut_seal","gmin","gnight", "g0")))


# Blackman, WTC4
wtc4gmin <- read.csv("data/wtc4_gmin_detached.csv") %>%
  rename(gmin = gmin_mmol_m2_s) %>%
  mutate(ch_temp_fac = as.factor(ch_temp))

wtc4gdark <- read.csv("data/wtc4_gnight.csv")




