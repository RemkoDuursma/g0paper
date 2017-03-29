
pacman::p_load(Hmisc, car, dplyr, tidyr, nlme, nlshelper, 
               forcats, tibble, magicaxis, 
               plantecophys, readxl, multcomp,
               reporttools)



if(!dir.exists("download"))dir.create("download")
if(!dir.exists("output"))dir.create("output")
source("R/functions.R")
source("R/figures.R")


# You must find TRY categorical traits yourself and place it in /data
# I cannot share it here (registration required)
if(!file.exists("data/trydb.rds")){
  tryfile <- "data/TRY_Categorical_Traits_Lookup_Table_2012_03_17_TestRelease.xlsx"
  if(!file.exists(tryfile))stop("Place the TRY look up table in data/.")
  trydb <- read_excel(tryfile) %>% 
    rename(species = AccSpeciesName)
  trydb <- trydb[,c("species","PhylogeneticGroup","PlantGrowthForm","LeafType","LeafPhenology")]
  saveRDS(trydb, "data/trydb.rds")
} else {
  trydb <- readRDS("data/trydb.rds")
}

# Additional growth form / phylogenetic group from googling missing species
gf_add <- read.csv("data/growthform_additional.csv", stringsAsFactors = FALSE)
trydb <- bind_rows(trydb, gf_add)



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
rieder <- read.csv("data/riederer2001_table1.csv", stringsAsFactors = FALSE) %>%
  filter(gmin < 10000) %>%
  mutate(gmin = 2 * 10^5 * 10^-3 * gmin / 41,
         source = "Riederer and Schreiber 2001",
         method = "gcut_isol")

# Kerstiens 1996
kerst <- read.csv("data/kerstiens1996_table1.csv", stringsAsFactors = FALSE) %>%
  filter(gmin > 0,
         species != "Molinia caerulea",
         method %in% c("A","C","D")) %>%
  mutate(method = dplyr::recode(method, A="gcut_isol", C="gcut_seal", D="gmin"),
         gmin = 2 * 10^5 * 10^-3 * gmin / 41,
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

# Phillips et al 2010
phillips <- read.csv("data/Phillips2010_gmin.csv", stringsAsFactors = FALSE) %>% 
  dplyr::select(species, gmin) %>%
  mutate(species = gsub("E. ","E.",species),
         species = gsub("E.","Eucalyptus ",species),
         source="Phillips et al. 2010") %>%
  group_by(species) %>%
  summarize(gmin = mean(gmin)) %>%
  mutate(species = as.character(species))

# Brodribb 2014
brodribb <- read.csv("data/brodribb2014_gmin_conifers.csv", stringsAsFactors = FALSE)
names(brodribb) <- c("family","species","gmin")
brodribb <- mutate(brodribb, gmin = 1000 * gmin,
                   source="Brodribb et al. 2014"
                   )

# gmin from various literature sources. 
gminrev <- read.csv("data/gmin_review_literature.csv",stringsAsFactors = FALSE) %>% 
  #mutate(method = "gmin") %>%
  #dplyr::select(gmin, species) %>% 
  mutate(gmin = 2 * 10^5 * 10^-3 * gmin / 41)

# Lopez, wet plants only
lopwet <- filter(lopez, treatment == "w") %>%
  mutate(source = "Rosana Lopez, unpublished") %>%
  group_by(species) %>%
  summarize(gmin = mean(gmin),
            source = first(source)) %>%
  mutate(species = as.character(species))

blackamb <- filter(wtc4gmin, growth_T == "amb") %>%
  dplyr::select(gmin) %>%
  mutate(species = "Eucalyptus parramattensis",
         source = "Chris Blackman, unpublished.") %>%
  group_by(species) %>%
  summarize(gmin = mean(gmin),
            source=first(source))

# Add brodribb, phillips, wet plants from Lopez, ambient plants from Blackman
gminrev <- bind_rows(gminrev, 
                     dplyr::select(brodribb, species, gmin, source),
                     phillips,
                     lopwet,
                     blackamb) %>%
  mutate(method = "gmin")

gminrev <- bind_rows(kerst, gminrev)

condreview_ave <- group_by(gminrev, species, method, source) %>%
  summarize(gmin = mean(gmin))

# Add plant growth form etc.
condreview_ave <- left_join(condreview_ave, trydb, by="species") %>%
  filter(!(PhylogeneticGroup %in% "Pteridophytes"))

condreview_ave$PlantGrowthForm <- as.factor(condreview_ave$PlantGrowthForm) %>%
  fct_collapse(shrub = c("herb/shrub","shrub","shrub/tree"),
               tree = "tree",
               herb = "herb",
               graminoid="graminoid")

condreview_ave$PhylogeneticGroup <- as.factor(condreview_ave$PhylogeneticGroup) %>%
  fct_collapse(Angiosperm = c("Angiosperm_Eudicotyl","Angiosperm_Magnoliid","Angiosperm_Monocotyl"))

condreview_ave$LeafPhenology <- as.factor(condreview_ave$LeafPhenology) %>%
  fct_collapse(evergreen = c("deciduous/evergreen","evergreen"))

condreview_ave <- ungroup(condreview_ave)

# Miner et al. 2016
miner <- read.csv("data/Miner_table1.csv")


# Lombardozzi et al 2017
lombar <- read.csv("data/lombardozzi_gnight.csv", stringsAsFactors = FALSE) %>%
  rename(gmin = gnight) %>%
  filter(gmin > 0) %>%
  mutate(method="gnight")


g0s <- data.frame(gmin=1000 * c(lin2015coef$g0, miner$g0),
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


gdfr <- bind_rows(ungroup(condreview_ave), lombar, g0s, minags) %>%
  mutate(method = factor(method, levels=c("gcut_isol","gcut_seal","gmin","gnight", "g0","gslowA")))







