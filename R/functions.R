to.pdf <- function(expr, filename, ..., verbose=TRUE) {
  if(!file.exists(dirname(filename)))
    dir.create(dirname(filename), recursive=TRUE)
  if ( verbose )
    cat(sprintf("Creating %s\n", filename))
  pdf(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}



fits_lin2015 <- function(lin2015a){
  
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
  
  return(p)
}


