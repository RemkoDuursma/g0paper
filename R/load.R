
if(!require("pacman", quietly=TRUE))install.packages("pacman")
suppressPackageStartupMessages(
  pacman::p_load(Hmisc, car, dplyr, tidyr, nlme, nlshelper, 
               forcats, tibble, magicaxis, 
               plantecophys, readxl, multcomp,
               reporttools,  tools, pander, knitr,
               doBy, stringi, rmarkdown, ggplot2)
)

if(!dir.exists("download"))dir.create("download")
if(!dir.exists("output"))dir.create("output")
source("R/functions.R")
source("R/figures.R")

# Lin et al. 2015
lin2015 <- download_and_read("https://ndownloader.figshare.com/files/1886204",
                             "download",
                             "lin2015.csv") %>%
  mutate(fitgroup = paste(Datacontrib, Species, sep="_"),
         Ci = CO2S - Photo/(Cond/1.6),
         CiCa = Ci / CO2S,
         BBopti = Photo/(CO2S*sqrt(VPD)))

# Filtered version, keep only groups with n>15, no bad PAR data, remove two groups.
lin2015a <- group_by(lin2015, fitgroup) %>%
  filter(n() > 15,
         PARin < 2500,
         PARin >= 0,
         Datacontrib != "Derek Eamus",  # sorry derek you are an outlier
         fitgroup != "Jon Bennie_Zornia")

# Lin 2015 fits: regression-based estimates of g0 and g1
lin2015coef <- fits_lin2015(lin2015a)

# Kerstiens 1996: literature compilation of various gmins.
kerst <- read.csv("data/kerstiens1996_table1.csv", stringsAsFactors = FALSE) %>%
  filter(gmin > 0,
         species != "Molinia caerulea",
         method %in% c("A","C")) %>%
  mutate(method = dplyr::recode(method, A="gcut_isol", C="gcut_seal"),
         gmin = 2 * 41 * 10^-5 * 10^3 * gmin,
         source = "Kerstiens 1996")
# Notes:
# A - astomatous cuticle, removed from leaf (permeance)
# B - detached leaves, stomatal side sealed
# C - whole leaves, stomatal side sealed (not clear when B or C, so all sealed leaves = C)
# D - mass loss of detached leaves (gmin)
# E1 - gnight / gdark
# E2 - minimum conductance measured during the day (vague!)
# E3 - presumed stomatal closure (v high CO2, ABA, drought, VPD) (vague!)

# Shuster et al. 2017 (JXB), only for gcut.
schuster2017 <- read.csv("data/schuster2017_tableS1.csv", stringsAsFactors = FALSE) %>%
  filter(parameter == "p") %>%
  rename(method = parameter, gmin = value) %>%
  mutate(method = "gcut_isol")


# Blackman, WTC4
wtc4gmin <- read.csv("data/wtc4_gmin_detached.csv") %>%
  rename(gmin = gmin_mmol_m2_s) %>%   # projected area
  mutate(ch_temp_fac = as.factor(ch_temp))

# WTC4 porometer data at night (should probably not trust)
wtc4gdark <- read.csv("data/wtc4_gnight.csv")

# Rosana Lopez Hakea data
lopez <- read.csv("data/lopez_gmin_hakea.csv") %>% 
  group_by(species, treatment) %>%
  dplyr::summarize(gmin = mean(gmin)) %>%
  mutate(gmin = 2*gmin)  # to convert - badly - to projected area

# some meta data for this dataset (leaf form, MAP)
lopmet <- read.csv("data/lopez_gmin_hakea_meta.csv")

# Wide format
lopw <- reshape(as.data.frame(lopez), direction="wide", idvar="species", timevar="treatment")

lopw <- dplyr::select(lopmet, species, leaf.form) %>% distinct() %>%
  right_join(lopw, by="species")

# Version 2; by population merged with lopmet (because MAP varies by population)
lopez2 <- read.csv("data/lopez_gmin_hakea.csv") %>% group_by(species, population, treatment) %>%
  dplyr::summarize(gmin = mean(gmin)) %>% 
  mutate(gmin = 2 * gmin) %>%
  inner_join(lopmet, by=c("species","population"))

# Miner et al. 2016
# Compilation of Ball-Berry parameters.
miner <- read.csv("data/Miner_table1.csv") %>%
  mutate(g0 = 1000*g0)


# Lombardozzi et al 2017
# Compilation of nighttime conductance.
lombar <- read.csv("data/lombardozzi_gnight.csv", stringsAsFactors = FALSE) %>%
  rename(gmin = gnight) %>%
  filter(gmin > 0, gmin < 1000, 
         method %in% c("ge","gas exchange","Li6400","IRGA")) %>%
  mutate(method="gnight")

# Combine Lin and Miner's estimate of g0 from 'regression' based.
g0s <- data.frame(gmin=c(lin2015coef$g0, miner$g0),
                  method="g0", stringsAsFactors = FALSE) %>%
  filter(!is.na(gmin),
         gmin > 0)


#------ gmindatabase
# Made with TRYDB and manual edits.
species_classes <- read.csv("data/Species_classifications.csv", stringsAsFactors = FALSE)

# Literature compilation of gmin data.
url <- paste0("https://raw.githubusercontent.com/RemkoDuursma/",
              "gmindatabase/master/combined/gmindatabase.csv")
gmindat <- download_and_read(url, "data") %>% 
  filter(gmin > 0) %>%
  group_by(species) %>%
  dplyr::summarize(gmin = mean(gmin),
            datasource = first(datasource)) %>%
  ungroup %>%
  left_join(species_classes, by="species")

# Crop data compilation (subset of the above, aggregated by genotype not species)
url <- paste0("https://raw.githubusercontent.com/RemkoDuursma/",
                            "gmindatabase/master/combined/cropgmindatabase.csv")
cropgmin <- download_and_read(url, "data")

# All data for gmin (i.e. raw compiled data, with all columns and rows for each study,
# not aggregated as the above two versions)
url <- paste0("https://raw.githubusercontent.com/",
              "RemkoDuursma/gmindatabase/master/combined/gminall.rds")
gminall <- download_and_read(url, "data")

# Read site climate data, prepared with the speciesmap R package.
# The binary intermediate file (climfile) is included in this repository,
# remove that file and run the code below to regenerate it.
# As the potential error shows, we need to install that package from github with:
# devtools::install_github("remkoduursma/speciesmap")
climfile <- "data/species_climate_wcpet.rds"
if(!file.exists(climfile)){
  
  r <- require(speciesmap, quietly=TRUE)
  if(!r){
    stop("Please install the speciesmap package like this (from R) :\n",
         "devtools::install_github('remkoduursma/speciesmap')\n",
         "if this fails, please check the README with this repository.\n",
         .call=FALSE
    )
  }
  
  dir.create("data/zomerpet", showWarnings = FALSE)
  dir.create("data/worldclim", showWarnings = FALSE)
  dir.create("data/ALAcache", showWarnings = FALSE)
  
  options(zomerpetpath="data/zomerpet", worldclimpath="data/worldclim")
  ALA4R::ala_config(cache_directory="data/ALAcache")
  
  sp <- unique(gmindat$species)
  nf <- sapply(strsplit(sp, " "), length)
  sp <- sp[nf == 2]
  sp <- sort(sp)
  
  wc <- climate_presence(sp, vars=c("tavg","prec","bio","pet"), database="both") %>%
    annualize_clim %>% aggregate_clim %>%
    mutate(species = as.character(species))
  
  saveRDS(wc, climfile)
} else {
  wc <- readRDS(climfile)
}

gmindat2 <- left_join(gmindat, wc, by="species")

# Add plant-growth form grouping variable (for one figure).
gmindat <- mutate(gmindat, group2 = case_when(Woody ~ "Woody",
                                              PlantGrowthForm == "graminoid" ~ "Graminoid",
                                              Crop ~ "Crop",
                                              TRUE ~ "Other non-woody"))


# Add taxonomic Family
# This binary file is included with this repository, 
# to regenerate it remove it and run this code again.
clsfile <- "data/cls_output.rds"
sp_dfr <- data.frame(species=unique(gmindat$species), stringsAsFactors = FALSE)
if(!file.exists(clsfile)){
  pacman::p_load(Taxonstand, taxize)
  cls <- classification(sp_dfr$species, db="ncbi")
  saveRDS(cls, clsfile)
} else {
  cls <- readRDS(clsfile)
}
cls <- as_dataframe_cls(cls)
rownames(cls) <- NULL
cls <- cls[!duplicated(cls$species),]

gmindat <- left_join(gmindat,  cls, by="species")

# gs at low PAR, gs at low A
# Both from updated Lin database. Extracted for convenience.
# See R_other/cond_low_a_low_par.R
gs_low_a <- read.csv("data/cond_low_a.csv") %>%
  dplyr::select(gmin = Cond.mean) %>%
  mutate(method = "gslowA",
         gmin = 1000*gmin)
gs_low_par <- read.csv("data/cond_low_par.csv") %>%
  dplyr::select(gmin = Cond) %>%
  mutate(method = "gslowPAR",
         gmin = 1000 * gmin)


# Dataframe with comparison methods.
gmindat_simple <- dplyr::select(gmindat, gmin) %>% 
  mutate(method = "gmin")

kerst_simple <- group_by(kerst, species, method) %>%
  dplyr::summarize(gmin = mean(gmin)) %>% 
  ungroup %>%
  dplyr::select(gmin, method)

gdfr <- bind_rows(gmindat_simple, kerst_simple, schuster2017, lombar, gs_low_a, gs_low_par) %>%
  filter(method != "gcut_seal") %>%  
  mutate(method = factor(method, levels=c("gcut_isol", "gmin", "gnight", "gslowPAR", "gslowA")))


# Herve Cochard's simulations with Sureau
planta <- read.csv("data/plant_a_sureau.csv", skip=1)
plantb <- read.csv("data/plant_b_sureau.csv", skip=1)


# Chris Blackman's example drying curves.
goodbadcurves <- read.csv("data/good and bad gmin.csv")



