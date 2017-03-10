
pacman::p_load(dplyr, tidyr, nlme, nlshelper, tibble)



if(!dir.exists("download"))dir.create("download")

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
         Datacontrib != "Derek Eamus",
         fitgroup != "Jon Bennie_Zornia")

