
kerst2a <- group_by(kerst2, species, method) %>%
  summarize(gmin=mean(gmin),
            PhylogeneticGroup = last(PhylogeneticGroup),
            PlantGrowthForm = last(PlantGrowthForm),
            LeafType = last(LeafType))

summ <- filter(kerst2a, 
                    method == "gmin",
                    !is.na(PlantGrowthForm)) %>%
  group_by(PlantGrowthForm) %>%
  summarize(gmin_mean = mean(gmin),
            n = n(),
            se = sd(gmin)/sqrt(n))


library(multcomp)

fit <- lm(gmin ~ PlantGrowthForm, data=kerst2a, subset=method=="gmin")
g <- glht(fit, linfct=mcp(PlantGrowthForm = "Tukey"))
cld(g)$mcletters$Letters


