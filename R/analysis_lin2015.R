

fits <- lmList(Cond ~ BBopti | fitgroup, data=lin2015a, na.action=na.omit)
p <- as.data.frame(coef(fits)) %>% 
  rownames_to_column(var = "fitgroup") %>%  
  rename(g0 = `(Intercept)`, g1 = BBopti) %>% 
  mutate(R2 = sapply(fits, function(x)summary(x)$adj.r.squared),
         rmse = sapply(fits, function(x)summary(x)$sigma),
         cor = sapply(fits, function(x)summary(x, correlation=T)$correlation[2,1]),
         g0_lci = sapply(fits, function(x)confint(x)[1,1]),
         g0_uci = sapply(fits, function(x)confint(x)[1,2])) %>% 
  left_join(summarize(lin2015a, 
                      sd_bbopti = sd(BBopti, na.rm=TRUE),
                      cv_bbopti = sd_bbopti/mean(BBopti, na.rm=TRUE),
                      sd_cond = sd(Cond, na.rm=TRUE),
                      Qrange = max(PARin) - min(PARin),
                      Drange = max(VPD) - min(VPD),
                      Trange = max(Tleaf) - min(Tleaf)), 
            by="fitgroup")

fits2 <- nlsList(Cond ~ g0 + 1.6*(1 + g1/sqrt(VPD))*(Photo / CO2S) | fitgroup, 
                 start=list(g0=0, g1=3),
                 data=lin2015a, na.action=na.omit)
p2 <- as.data.frame(coef(fits2)) %>% 
  rename(g0f = g0, g1f = g1) %>%
  mutate(rmsef = sapply(fits2, function(x)summary(x)$sigma))
p <- cbind(p, p2)

fits3 <- nlsList(Cond ~ 1.6*(1 + g1/sqrt(VPD))*(Photo / CO2S) | fitgroup, 
                 start=list(g1=3),
                 data=lin2015a, na.action=na.omit)
p3 <- as.data.frame(coef(fits3)) %>% 
  rename(g1n = g1) %>%
  mutate(rmsen = sapply(fits3, function(x)summary(x)$sigma))
p <- cbind(p, p3)



hist(p$g0, breaks=30)




par(mfrow=c(1,2))
l <- loess(g0 ~ R2, data=p, span=0.8)
plot_loess(l, pch=16, col="dimgrey")
abline(h=0, lty=3)

with(p, plot(g1, g0))
add_regres_line(lm(g0 ~ g1, data=p))
abline(h=0, lty=3)



with(p, {
  plot(R2, g0, pch=19,
       panel.first=segments(x0=R2, x1=R2, y0=g0_lci, y1=g0_uci, col="dimgrey")
  )
})
abline(h=0)
l <- loess(g0 ~ R2, data=p, span=0.8)
plot_loess(l, add=TRUE, band=FALSE, lwd=2, lty=3)




minags <- group_by(lin2015, fitgroup) %>%
  summarize(
               Amin = min(Photo, na.rm=TRUE),
               gsmin = min(Cond, na.rm=TRUE),
               Qrange = max(PARin) - min(PARin),
               Pathway = unique(Pathway)[1]
               ) %>%
  filter(Amin < 5)

with(minags, plot(Amin, gsmin, pch=19, col=ifelse(Pathway == "C3", "blue","red")))






