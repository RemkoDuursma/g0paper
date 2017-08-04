
pacman::p_load(Hmisc, car, dplyr, tidyr, nlme, nlshelper, 
               forcats, tibble, magicaxis, 
               plantecophys, readxl, multcomp,
               reporttools, Taxonstand, taxize, speciesmap)



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


# Blackman, WTC4
wtc4gmin <- read.csv("data/wtc4_gmin_detached.csv") %>%
  rename(gmin = gmin_mmol_m2_s) %>%   # projected area
  mutate(ch_temp_fac = as.factor(ch_temp))

# WTC4 porometer data at night (should probably not trust)
wtc4gdark <- read.csv("data/wtc4_gnight.csv")

# Rosana Lopez Hakea data
lopez <- read.csv("data/lopez_gmin_hakea.csv") %>% group_by(species, treatment) %>%
  summarize(gmin = mean(gmin)) %>%
  mutate(gmin = 2*gmin)  # to convert - badly - to projected area

# some meta data for this dataset (leaf form, MAP)
lopmet <- read.csv("data/lopez_gmin_hakea_meta.csv")

# Wide format
lopw <- reshape(as.data.frame(lopez), direction="wide", idvar="species", timevar="treatment")

lopw <- dplyr::select(lopmet, species, leaf.form) %>% distinct() %>%
  right_join(lopw, by="species")

# Version 2; by population merged with lopmet (because MAP varies by population)
lopez2 <- read.csv("data/lopez_gmin_hakea.csv") %>% group_by(species, population, treatment) %>%
  summarize(gmin = mean(gmin)) %>% 
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
  filter(gmin > 0) %>%
  mutate(method="gnight")

# Combine Lin and Miner's estimate of g0 from 'regression' based.
g0s <- data.frame(gmin=c(lin2015coef$g0, miner$g0),
                  method="g0", stringsAsFactors = FALSE) %>%
  filter(!is.na(gmin),
         gmin > 0)


# From Lin2015, data where A < threshold.
minags <- group_by(lin2015, fitgroup) %>%
  summarize(
    Amin = min(Photo, na.rm=TRUE),
    gmin = 1000 * min(Cond, na.rm=TRUE),
    Qrange = max(PARin) - min(PARin),
    Pathway = unique(Pathway)[1]
  ) %>%
  filter(Amin < 2) %>%
  mutate(method="gslowA") %>%
  dplyr::select(gmin, method)


#------ gmindatabase
# Made with TRYDB and manual edits.
species_classes <- read.csv("data/Species_classifications.csv", stringsAsFactors = FALSE)

# Data from gmindatabase
# gmindat <- read.csv("https://raw.githubusercontent.com/RemkoDuursma/gmindatabase/master/combined/gmindatabase.csv",
#                     stringsAsFactors = FALSE) 
gmindat <- read.csv("c:/repos/gmindatabase/combined/gmindatabase.csv",
                    stringsAsFactors = FALSE) %>%
  filter(gmin > 0) %>%
  group_by(species) %>%
  summarize(gmin = mean(gmin),
            source = first(source)) %>%
  ungroup %>%
  left_join(species_classes, by="species")

#
climfile <- "data/species_climate_wcpet.rds"
if(!file.exists(climfile)){
  options(zomerpetpath="c:/data/zomerpet", worldclimpath="c:/data/worldclim")
  ALA4R::ala_config(cache_directory="c:/data/ALAcache")
  
  sp <- unique(gmindat$species)
  nf <- sapply(strsplit(sp, " "), length)
  sp <- sp[nf == 2]
  sp <- sort(sp)
  
  wc <- climate_presence(sp, vars=c("tavg","prec","bio","pet"), database="both")
  
  saveRDS(wc, climfile)
} else {
  wc <- readRDS(climfile)
}

wc2 <- annualize_clim(wc) %>% aggregate_clim %>%
  mutate(species = as.character(species))

gmindat2 <- left_join(gmindat, wc2, by="species")

# New grouping variable.
gmindat$group2 <- NA
gmindat$group2[gmindat$Woody] <- "Woody"
gmindat$group2[gmindat$PlantGrowthForm == "graminoid"] <- "Graminoid"
gmindat$group2[gmindat$Crop] <- "Crop"
gmindat$group2[is.na(gmindat$group2)] <- "Other non-woody"


# Add Family
clsfile <- "data/cls_output.rds"
sp <- unique(gmindat$species)
if(!file.exists(clsfile)){
  cls <- classification(sp, db="ncbi")
  saveRDS(cls, clsfile)
  } else {
    cls <- readRDS(clsfile)
  }

ord <- sapply(cls, function(x)if(!all(is.na(x)))x$name[x$rank == "order"] else NA) %>% unname
fam <- sapply(cls, function(x)if(!all(is.na(x)))x$name[x$rank == "family"] else NA) %>% unname
gmindat <- left_join(gmindat, 
                     data.frame(species=sp, Order=ord, Family=fam, stringsAsFactors = FALSE),
                     by="species")


# Dataframe with comparison methods.
gmindat_simple <- dplyr::select(gmindat, gmin) %>% 
  mutate(method = "gmin")

kerst_simple <- group_by(kerst, species, method) %>%
  summarize(gmin = mean(gmin)) %>% 
  ungroup %>%
  dplyr::select(gmin, method)

gdfr <- bind_rows(gmindat_simple, kerst_simple, lombar, g0s, minags) %>%
  mutate(method = factor(method, levels=c("gcut_isol","gcut_seal","gmin","gnight", "g0","gslowA")))







