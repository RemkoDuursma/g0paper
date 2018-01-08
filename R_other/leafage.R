

ii <- which(sapply(gminall, function(x)"leafage" %in% names(x)))
datsub <- gminall[ii]

sel <- function(x){
  x <- x[,c("species","leafage","gmin","datasource")]
  x$leafage <- as.character(x$leafage)
  x
}
datsub <- lapply(datsub, sel)

dat <- bind_rows(datsub) %>%
  group_by(species, leafage, datasource) %>%
  summarize(gmin = mean(gmin, na.rm=TRUE)) %>%
  ungroup


dat$leafage[dat$leafage %in% c("current year","current")] <- "0"
dat$leafage[dat$leafage %in% c("one year old","one-year-old")] <- "1"
dat$leafage[dat$leafage == "C"] <- "0"
dat$leafage[dat$leafage == "C+1"] <- "1"
dat$leafage[dat$leafage == "C+2"] <- "2"
dat$leafage[dat$leafage == "C+3"] <- "3"

dat$leafage <- factor(dat$leafage, levels=c("young","old","0","1","2","3","4","5","10"))


library(ggplot2)


ggplot(dat, aes(x=leafage, y=gmin)) +
  geom_line(aes(group=species)) + 
  geom_point() + 
  geom_hline(yintercept=0, alpha=0) + 
  facet_wrap(~datasource, drop=TRUE, scales= "free") +
  theme_bw() +
  scale_colour_manual(values = rainbow(12)) +
  labs(x="", y=expression(g[min]~~(mmol~m^-2~s^-1)))


