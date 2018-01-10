
library(lme4)
library(lmerTest)

fit1 <- lmer(log(gmin) ~ growth_T + factor(ch_temp) + (1|chamber), data=wtc4gmin)
anova(fit1, test = "F")
car::Anova(fit1)

group_by(wtc4gmin, growth_T) %>%
  summarize(gmin = mean(gmin))



group_by(gdfr, method) %>%
  summarize(gmin = exp(mean(log(gmin))))

table(gdfr$method)
