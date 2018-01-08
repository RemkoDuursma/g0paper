
choat <- read_excel("data/xylem functional traits database Master May 2014.xls") %>%
  mutate(species = paste(Genus, Species)) %>%
  rename(P50 =`P50 (MPa)`) %>%
  dplyr::select(species, P50) %>%
  filter(!is.na(P50)) %>%
  group_by(species) %>%
  summarize(P50 = mean(P50))



gmindat2 <- inner_join(gmindat, choat, by="species")
with(gmindat2, plot(P50, gmin, 
                    pch=c(1,19)[as.factor(PhylogeneticGroup)]))

with(subset(gmindat2, datasource != "brodribb2014"),
     plot(P50, gmin, 
                    pch=c(1,19)[as.factor(PhylogeneticGroup)]))


hist(choat$P50)

hist(log(-choat$P50))
hist(log(martin$P50m))


# Sample Pgs90 from mean function with P50
# plus SD of residuals.
# See gmindat_vs_martin
sample_pgs90 <- function(P50){
  3.39*(1 - exp(-0.282*P50)) +
    rnorm(1, mean = 0, sd=0.568)
}

# Sample p50 from truncated log-normal
sample_p50 <- function(){
  draw <- 20
  while(draw > 15 | draw < 2){
    draw <- rlnorm(1, meanlog = 1.67, sdlog = 0.51)
  }
return(draw)
}

sample_p50_pgs90 <- function(n){

  p50s <- replicate(n, sample_p50())
  
  pgs90 <- sapply(p50s, function(x)sample_pgs90(x))
  
return(data.frame(P50=-p50s, Pgs90=-pgs90))
}


sample_gmin <- function(){
  draw <- 30
  while(draw > 25 | draw < 0.2){
    draw <- rlnorm(1, meanlog = 1.29, sdlog = 0.85)
  }
  return(draw)
}

#check

with(martin, plot(P50, Pgs90, ylim=c(0,-5)))
with(sample_p50_pgs90(nrow(martin)), 
     points(P50, Pgs90, pch=19))




gmintree <- subset(gmindat, PlantGrowthForm == "tree")
mean(log(gmintree$gmin))
sd(log(gmintree$gmin))

