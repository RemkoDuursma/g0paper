
pacman::p_load(car, dplyr, tidyr, nlme, nlshelper, 
               forcats, tibble, magicaxis, 
               plantecophys, readxl)



if(!dir.exists("download"))dir.create("download")
if(!dir.exists("output"))dir.create("output")
source("R/functions.R")
source("R/figures.R")


# You must find TRY categorical traits yourself and place it in /data
# I cannot share it here (registration required)
if(!file.exists("data/trydb.rds")){
  tryfile <- "data/TRY_Categorical_Traits_Lookup_Table_2012_03_17_TestRelease.xlsx"
  if(!file.exists(tryfile))stop("Place the TRY look up table in data/.")
  trydb <- read_excel(tryfile)
  trydb <- trydb[,c("AccSpeciesName","PhylogeneticGroup","PlantGrowthForm","LeafType","LeafPhenology")]
  saveRDS(trydb, "data/trydb.rds")
} else {
  trydb <- readRDS("data/trydb.rds")
}



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


# Riederer and Schreiber 2001
rieder <- read.csv("data/riederer2001_table1.csv") %>%
  mutate(gmin = 10^5 * 10^-3 * gmin / 41)

# Kerstiens 1996
kerst <- read.csv("data/kerstiens1996_table1.csv", stringsAsFactors = FALSE) %>%
  filter(gmin > 0,
         method %in% c("A","C","D")) %>%
  mutate(method = dplyr::recode(method, A="gcut_isol", C="gcut_seal", D="gmin"),
         gmin = 2 * 10^5 * 10^-3 * gmin / 41)   # convert to mmol m-2 s-1 assume 41 mol m-3
                # 2 because gmin was expressed as per total leaf area.

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


# Kerstiens dataset with growth form, phenology, etc.
kerst2 <- inner_join(kerst, trydb, by=c("species" = "AccSpeciesName")) %>%
  filter(PhylogeneticGroup != "Pteridophytes")

kerst2$PlantGrowthForm <- as.factor(kerst2$PlantGrowthForm) %>%
  fct_collapse(shrub = c("herb/shrub","shrub","shrub/tree"),
               tree = "tree",
               herb = "herb",
               graminoid="graminoid")

kerst2$PhylogeneticGroup <- as.factor(kerst2$PhylogeneticGroup) %>%
  fct_collapse(Angiosperm = c("Angiosperm_Eudicotyl","Angiosperm_Magnoliid","Angiosperm_Monocotyl"))


kerst2$LeafPhenology <- as.factor(kerst2$LeafPhenology) %>%
  fct_collapse(evergreen = c("deciduous/evergreen","evergreen"))


# Blackman, WTC4
wtc4gmin <- read.csv("data/wtc4_gmin_detached.csv") %>%
  rename(gmin = gmin_mmol_m2_s) %>%
  mutate(ch_temp_fac = as.factor(ch_temp))

wtc4gdark <- read.csv("data/wtc4_gnight.csv")


# From Lin2015, data where A < threshold.
minags <- group_by(lin2015, fitgroup) %>%
  summarize(
    Amin = min(Photo, na.rm=TRUE),
    gsmin = min(Cond, na.rm=TRUE),
    Qrange = max(PARin) - min(PARin),
    Pathway = unique(Pathway)[1]
  ) %>%
  filter(Amin < 2)



