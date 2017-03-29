
library(stringr)

phil <- as.data.frame(read_excel("data/Cuticular Conductance 6-1-10 (RA).xls", skip=7))

names(phil)[c(1:2, ncol(phil))] <- c("species", "leafage", "leafarea")
phil <- dplyr::select(phil, species, leafage, leafarea, everything())

# time information into columns; useful for reshape
names(phil)[4:ncol(phil)] <- paste0("time.", phil[1,4:ncol(phil)])


phil <- slice(phil, -c(1:2)) %>%
  filter(!is.na(species))
  

phill <- reshape(as.data.frame(phil), varying=4:ncol(phil), direction="long", timevar="Time") %>%
  rename(weight=time) %>% 
  mutate(species = str_trim(species),
         leafage = str_trim(leafage)) %>%
  arrange(id, Time)


#for(i in 1:30)with(subset(phill, id == i), plot(Time, weight, main=i))

badcurves <- c(3,7,23,25,28,29)

phillsub <- filter(phill, 
                   !id %in% badcurves,
                   Time < 1400)



library(segmented)
x <- subset(phillsub, id ==1)
fit <- lm(weight ~ Time, data=x)
seg <- segmented(fit, seg.Z = ~Time)


for(i in unique(phillsub$id)){
  x <- subset(phillsub, id == i)
  fit <- lm(weight ~ Time, data=x)
  seg <- segmented(fit, seg.Z = ~Time)
  
  with(x, plot(Time, weight, main=i))
  plot(seg, add=TRUE)
}


phils <- split(phillsub, phillsub$id)
fits <- lapply(phils, function(x){
  fit <- lm(weight ~ Time, data=x)
  return(segmented(fit, seg.Z = ~Time))
})

p <- sapply(fits, coef)

# coefficients are intercept, slope section 1, slope section 2 as difference from section 1
# don't know what coefficient 4 is but always zero

# emin in g min-1
emin <- apply(p, 2, function(x)x[2] + x[3])

emindf <- data.frame(id = as.numeric(names(emin)), emin=-as.vector(emin))

m <- dplyr::select(phillsub, id, species, leafage, leafarea) %>% distinct()

emindf <- left_join(emindf, m, by="id") %>%
  mutate(eleaf = (1/18) * (10^3 * emin / 60) / (leafarea * 10^-4))


# RH was 50.7%, Tair 25C though both fluctuated quite a bit!!
VPD <- plantecophys::RHtoVPD(50.7, 25)

emindf$gmin <- emindf$eleaf / (VPD / 101)

# export
write.csv(emindf[,c("species","leafage","eleaf","leafarea","gmin")],
          "data/Phillips2010_gmin.csv", row.names=FALSE)

write.csv(phill, "data/Phillips2010_gminrawdata.csv", row.names=FALSE)




