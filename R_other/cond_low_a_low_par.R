
# This works only in the repository bitbucket.org/bmedlyn/leafgasexchange_ongoing

source("make_database.R")

library(doBy)
library(dplyr)

dat <- filter(alldata, opt %in% c("opt","low PAR","lower canopy"),
              PARin > 0,
              PARin < 40) %>%
  mutate(dataset_fitgroup = paste(Datacontrib, Location, fitgroup))


cond_low_par <- summaryBy(Cond ~ dataset_fitgroup, FUN=mean, data=dat, keep.names=TRUE)

dat <- filter(alldata, PARin > 40, Photo < 1, opt == "opt") %>%
  mutate(dataset_fitgroup = paste(Datacontrib, Location, fitgroup))
cond_low_a <- summaryBy(Cond ~ dataset_fitgroup, FUN=mean, data=dat)


write.csv(cond_low_par, "cond_low_par.csv", row.names = FALSE)
write.csv(cond_low_a, "cond_low_a.csv", row.names = FALSE)
