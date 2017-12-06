# useless idea; only 33 species overlap   


#https://zenodo.org/record/854700#.Whyx7kqWbIV

martin <- read_excel("data/DataBase.xlsx", sheet=4) %>%
  rename(species = Species.binomial)
names(martin)[6] <- "Pgs90"

martin_sub <- dplyr::select(martin, species, P50, Pgs90)

gmindat2 <- inner_join(gmindat, martin_sub, by="species")

with(gmindat2, plot(Pgs90, gmin))



martin <- read.csv("c:/repos/nlshelper/R_other/martin.csv")
library(nlshelper)
martin$P50m <- -martin$P50
martin$Pgs90m <- -martin$Pgs90
fit1 <- nls(Pgs90m ~ ymax*(1-exp(-k*P50m)), 
            start=list(ymax=2, k=1), data=martin)
plot_nls(fit1)







martin <- read_excel("data/DataBase.xlsx", sheet=4) %>%
  rename(species = Species.binomial)
names(martin)[6] <- "Pgs90"

martin <- mutate(martin, 
                 P50m = -P50,
                 Pgs90m = -Pgs90) %>%
  filter(!is.na(Pgs90))


fit2 <- nls(Pgs90m ~ SSasympOff(P50m, Asym, lrc, c0), data=martin)
plot_nls(fit2, xlim=c(0,20), ylim=c(0,5))


library(quantreg)

fit_nl1 <- lapply(c(0.1,0.5,0.9),
                  function(x)nlrq(Pgs90m ~ SSasympOff(P50m, Asym, lrc, c0), 
                       data=martin,
                       tau=x))
                

plot_nls(fit_nl1[[2]], xlim=c(0,20), ylim=c(0,5))
plot_nls(fit_nl1[[1]], add=T, lty=5)
plot_nls(fit_nl1[[3]], add=T, lty=5)

