

figure_gmin_review <- function(gdfr){
  
  par(mar=c(3,5,1,0.5), cex.lab=1.2)
  fit <- lm(log10(gmin) ~ method-1, data=gdfr)
  cis <- confint(fit)
  quans <- sapply(split(gdfr, gdfr$method), function(x)quantile(log10(x$gmin), probs=c(0.05, 0.95)))
  gmins <- with(gdfr, tapply(log10(gmin), method, mean, na.rm=TRUE))
  
  plot(1:5, gmins, pch=19, cex=1.2, axes=FALSE, xlab="",
       panel.first=segments(x0=1:5, x1=1:5, y0=quans[1,], y1=quans[2,], col="grey"),
       ylim=c(0,2),
       xlim=c(0.5,5.5),
       ylab=expression(Conductance~~(mmol~m^-2~s^-1)))
  axis(1, at=1:5, labels=c(expression(g["cut,isol"]), 
                           expression(g["min,seal"]), 
                           expression(g[min]), 
                           expression(g[dark]),
                           expression(g[0])))
  magaxis(side=2, unlog=2)
  box()
  arrows(x0=1:5, x1=1:5, y0=cis[,1], y1=cis[,2], angle=90, length=0.1, code=3)
  
  
}

figure_gmin_review_2 <- function(gdfr, minags){
  par(mar=c(3,5,1,0.5), cex.lab=1.2)
  fit <- lm(log10(gmin) ~ method-1, data=gdfr)
  cis <- confint(fit)
  quans <- sapply(split(gdfr, gdfr$method), function(x)quantile(log10(x$gmin), probs=c(0.05, 0.95)))
  gmins <- with(gdfr, tapply(log10(gmin), method, mean, na.rm=TRUE))
  
  plot(1:5, 10^gmins, pch=19, cex=1.2, xlab="",
       axes=FALSE,
       panel.first=segments(x0=1:5, x1=1:5, y0=10^quans[1,], y1=10^quans[2,], col="grey"),
       ylim=c(0,50),
       xlim=c(0.5,6.5),
       ylab=expression(Conductance~~(mmol~m^-2~s^-1)))
  axis(1, at=1:6, labels=c(expression(g["cut,isol"]),
                           expression(g["min,seal"]),
                           expression(g[min]),
                           expression(g[dark]),
                           expression(g[0]),
                           expression(g["s,low A"])))
  axis(2)
  box()
  arrows(x0=1:5, x1=1:5, y0=10^cis[,1], y1=10^cis[,2], angle=90, length=0.1, code=3)
  
  quang <- quantile(minags$gsmin, probs=c(0.05, 0.95))
  points(6, 1000*mean(minags$gsmin), pch=19, cex=1.2,
         panel.first=segments(x0=6, x1=6, y0=1000*quang[1], y1=1000*quang[2], col="grey"))
  cig <- 1000*(mean(minags$gsmin) + c(-2,2)*sd(minags$gsmin)/sqrt(nrow(minags)))
  arrows(x0=6, x1=6, y0=cig[1], y1=cig[2], angle=90, length=0.1, code=3)
}


figure_gmin_review_3 <- function(kerst2){
  
  plot_kerst2 <- function(yvar = "PhylogeneticGroup", meth = "gcut_isol", data=kerst2, ...){
    
    data <- data[data$method == meth,]
    data$yvar <- as.factor(data[,yvar])
    
    fit <- lm(log10(gmin) ~ yvar - 1, data=data)
    ci <- 10^confint(fit)
    mn <- 10^coef(fit)
    
    g <- glht(fit, linfct=mcp(yvar = "Tukey"))
    lets <- cld(g)$mcletters$Letters
    
    n <- length(mn)
    plot(1:n, mn, pch=19, cex=1.1, ylim=c(0, max(ci)+0.1*max(ci)), axes=FALSE,
         xlim=c(0.7, n+0.3), ...)
    axis(1, at=1:n, labels=capitalize(levels(data$yvar)))
    arrows(x0=1:n, x1=1:n, y0=ci[,1], y1=ci[,2], angle=90, length=0.05, code=3)
    text(x=1:n, y=ci[,2], lets, pos=3)
    axis(2)
    box()
  }
  
  # plot_kerst2("PlantGrowthForm", "gmin")
  # plot_kerst2("PhylogeneticGroup", "gmin", subset(kerst2, PlantGrowthForm == "tree"))
  # #plot_kerst2("PhylogeneticGroup", "gmin")
  
  # Combine angio/gymno with growth form
  kerst2$group <- as.character(kerst2$PlantGrowthForm)
  kerst2$group[which(kerst2$PlantGrowthForm == "tree")] <- 
    as.character(kerst2$PhylogeneticGroup)[which(kerst2$PlantGrowthForm == "tree")]
  kerst2$group <- factor(kerst2$group, levels=c("graminoid","herb","shrub","Angiosperm","Gymnosperm"))
  
  par(mar=c(6.5,4,1,1), mgp=c(2.5,0.5,0), tcl=-0.2, 
      las=3, yaxs="i",
      cex.lab=1.2, cex.axis=0.8)
  plot_kerst2("group", "gmin",
              xlab="", ylab=expression(g[min]~~(mmol~m^-2~s^-1)))
  #mtext(side=1, text="Tree", at=4.5, line=2, cex=1)
  
}




figure_g0g1_cor <- function(lin2015){
  
  x <- subset(lin2015, fitgroup == "Nicolas Martin-StPaul_Quercus ilex_StPaul_Puechabon")
  
  # for bootCase to work
  assign("x",x,envir=.GlobalEnv)
  
  fit <- lm(Cond ~ BBopti, data=x)
  
  el <- confidenceEllipse(fit, draw=FALSE)
  set.seed(123)
  
  b <- bootCase(fit)
  rm(x, envir=.GlobalEnv)
  
  par(mfrow=c(1,2), mar=c(4,4,1,1), mgp=c(2.5,0.5,0), tcl=0.2, las=1,
      cex.lab=1.1, cex.axis=0.9)
  with(x, plot(BBopti, Cond, ylim=c(0, 0.15),
               xlim=c(0,0.04),
               panel.first=add_regres_line(fit),
               xlab=expression(A/(C[a]*sqrt(D))),
               ylab=expression(g[s]~(mol~m^-2~s^-1)),
               pch=16, cex=0.8, col="dimgrey"))
  box()
  
  plot(el, type='l', lty=3, 
       xlab=expression(g[0]~~(mol~m^-2~s^-1)), 
       ylab=expression(g[1]~~(kPa^-0.5)),
       ylim=c(3,5))
  points(b, pch=16, cex=0.3, col="dimgrey")
  points(coef(fit)[1],coef(fit)[2], pch=19)
  
}




figure_R2g0 <- function(lin2015coef, miner){
  
  par(mfrow=c(1,2), mar=c(5,5,1,1), cex.axis=0.9)
  with(lin2015coef, {
    plot(R2, g0, pch=16,
         xlab=expression(R^2),
         ylab=expression(g[0]~~(mol~m^-2~s^-1)),
         ylim=c(-0.1, 0.3),
         xlim=c(0,1),
         panel.first=segments(x0=R2, x1=R2, y0=g0_lci, y1=g0_uci, col="dimgrey")
    )
  })
  abline(h=0)
  l <- loess(g0 ~ R2, data=lin2015coef, span=0.8)
  plot_loess(l, add=TRUE, band=TRUE, lwd=2, col="darkgrey")
  legend("topright", "Lin et al. 2015", bty='n', cex=0.8)
  
  with(miner,plot(R2, g0, pch=16, ylim=c(-0.1, 0.3), xlim=c(0,1),
                  xlab=expression(R^2),
                  ylab=expression(g[0]~~(mol~m^-2~s^-1))))
  abline(h=0)
  l <- loess(g0 ~ R2, data=miner, span=0.8)
  plot_loess(l, add=TRUE, band=TRUE, lwd=2, col="darkgrey")
  legend("topright", "Miner et al. 2017", bty='n', cex=0.8)
  
}


figure_amings <- function(lin2015){
  
  par(mar=c(5,5,1,1), cex.axis=0.9, cex.lab=1.1)

  with(minags, plot(Amin, gsmin, pch=16,
                    ylim=c(0,0.15),
                    xlab=expression(min~A~~(mu*mol~m^-2~s^-1)),
                    ylab=expression(min~g[s]~(mol~m^-2~s^-1)),
                    panel.first=add_regres_line(lm(gsmin~Amin))))
  
}


figure_sim <- function(){
  
  par(mfrow=c(2,2), mar=c(4,4,1,1), mgp=c(2.5,0.5,0), tcl=0.2, las=1,
      cex.axis=0.9, cex.lab=1.1)
  
  g1 <- 4
  g0 <- 0.03
  
  lty0 <- 1
  lty1 <- 5
  col0 <- "black"
  col1 <- "black"
  
  r0 <- Photosyn(PPFD=seq(20,1000, length=101),
                 g1=g1, g0=0, Vcmax=70, Jmax=140,
                 VPD=2, Tleaf=20, Ca=400)
  r1 <- Photosyn(PPFD=seq(20,1000, length=101),
                 g1=g1, g0=g0, Vcmax=70, Jmax=140,
                 VPD=2, Tleaf=20, Ca=400)
  
  with(r0, plot(PPFD, ALEAF/GS, type='l', ylim=c(0,80), col=col0, lty=lty0,
                xlab=expression(PPFD~~(mu*mol~m^-2~s^-1)),
                ylab=expression(A[n]/g[s]~~(mu*mol~mol^-1))
                ))
  with(r1, lines(PPFD, ALEAF/GS, col=col1, lty=lty1))
  
  legend("bottomright", c(expression(g[0] == 0), expression(g[0] == 0.03)),
         lty=c(lty0, lty1), bty='n')
  
  with(r0, plot(PPFD, Ci, type='l',
                xlab=expression(PPFD~~(mu*mol~m^-2~s^-1)),
                ylab=expression(C[i]~~(mu*mol~mol^-1)),
                col=col0, lty=lty0,
                ylim=c(250,400)))
  with(r1, lines(PPFD, Ci, col=col1, lty=lty1))
  
  tleafs <- seq(10,45,length=101)
  vpds <- 0.000605 * tleafs^2.39
  r0 <- Photosyn(PPFD=1500,
                 g1=g1, g0=0, Vcmax=70, Jmax=140,
                 VPD=vpds, Tleaf=tleafs, Ca=400)
  r1 <- Photosyn(PPFD=1500,
                 g1=g1, g0=g0, Vcmax=70, Jmax=140,
                 VPD=vpds, Tleaf=tleafs, Ca=400)
  
  with(r0, plot(Tleaf, ELEAF, type='l', 
                xlab=expression(T[leaf]~~(degree*C)),
                ylab=expression(E[L]~~(mmol~m^-2~s^-1)),
                col=col0, lty=lty0,
                ylim=c(0,10)))
  with(r1, lines(Tleaf, ELEAF, col=col1, lty=lty1))
  
  with(r0, plot(VPD, Ci, type='l', 
                xlab="VPD (kPa)",
                ylab=expression(C[i]~~(mu*mol~mol^-1)),
                col=col0, lty=lty0,
                ylim=c(0,400)))
  with(r1, lines(VPD, Ci, col=col1, lty=lty1))
  
}





figure_wtc4_gmin <- function(wtc4gmin, wtc4gdark){
  
  l <- layout(matrix(c(1,2), ncol=2), widths=c(2,1))
  
  par(mar=c(4,4,1,1), mgp=c(2.5, 0.5, 0), tcl=0.1, cex.axis=0.9)
  wtc4gmina <- group_by(wtc4gmin, chamber, ch_temp) %>%
    summarize(gmin = mean(gmin),
              growth_T = unique(growth_T)) %>%
    group_by(ch_temp, growth_T) %>%
    summarize(gmin_mean = mean(gmin),
              n = n(),
              se = sd(gmin)/sqrt(n),
              gmin_lcl = gmin_mean - qt(0.975, n-1)*se,
              gmin_ucl = gmin_mean + qt(0.975, n-1)*se
    )
  
  palette(c("blue2", "red2"))
  
  with(wtc4gmina, {
    plot(ch_temp, gmin_mean, ylim=c(0,12), pch=19, col=growth_T, cex=1.2,
         xlab=expression(Measurement~T~~(degree*C)),
         ylab=expression(g[min]~~(mmol~m^-2~s^-1)),
         xlim=c(17,28), axes=FALSE)
    arrows(x0=ch_temp, x1=ch_temp, y0=gmin_mean - se, y1=gmin_mean + se, angle=90, code=3, length=0.07, col=growth_T)
  })
  axis(2)
  axis(1, at=seq(17.5, 27.5, by=2.5))
  box()
  legend("topleft", c("Ambient",expression(Ambient~+~3~degree*C)), pch=19, pt.cex=1.1, col=palette(), 
         title="Growth T", bty='n')
  
  
  # library(lme4)
  # fit <- lmer(gmin ~ ch_temp + growth_T -1 + (1|chamber), data=wtc4gmin)
  # Anova(fit, test="F")
  
  wtc4gdarka <- group_by(wtc4gdark, chamber, surface) %>%
    summarize(gdark = mean(g_night),
              growth_T = unique(trt)) %>%
    group_by(surface, growth_T) %>%
    summarize(gdark_mean = mean(gdark),
              n = n(),
              se = sd(gdark)/sqrt(n)
    )
  
  xv <- c(0.88, 1.12, 1.88, 2.12)
  with(wtc4gdarka, {
    plot(xv, gdark_mean, ylim=c(0,300), axes=FALSE, pch=19, cex=1.2, col=growth_T,
         ylab=expression(g[dark]~~(mmol~m^-2~s^-1)),
         xlim=c(0.65, 2.45),
         xlab="")
    arrows(x0=xv, x1=xv, y0=gdark_mean - se, y1=gdark_mean + se, angle=90, code=3, length=0.07, col=growth_T)
  })
  axis(2)
  par(mgp=c(2.5,1.5, 0))
  axis(1, at=c(mean(xv[1:2]), mean(xv[3:4])), labels=c("Lower\nsurface", "Upper\nsurface"))
  box()
  
}



figure_wtc4_gmin_2 <- function(wtc4gmin){
  
  #l <- layout(matrix(c(1,2), ncol=2), widths=c(2,1))
  
  par(mar=c(4,4,1,1), mgp=c(2.5, 0.5, 0), tcl=0.2, cex.axis=0.9, 
      cex.lab=1.2, yaxs="i")
  wtc4gmina <- group_by(wtc4gmin, chamber, ch_temp) %>%
    summarize(gmin = mean(gmin),
              growth_T = unique(growth_T)) %>%
    group_by(ch_temp, growth_T) %>%
    summarize(gmin_mean = mean(gmin),
              n = n(),
              se = sd(gmin)/sqrt(n),
              gmin_lcl = gmin_mean - qt(0.975, n-1)*se,
              gmin_ucl = gmin_mean + qt(0.975, n-1)*se
    )
  
  palette(c("blue2", "red2"))
  
  with(wtc4gmina, {
    plot(ch_temp, gmin_mean, ylim=c(0,12), pch=19, col=growth_T, cex=1.2,
         xlab=expression(Measurement~T~~(degree*C)),
         ylab=expression(g[min]~~(mmol~m^-2~s^-1)),
         xlim=c(17,28), axes=FALSE)
    arrows(x0=ch_temp, x1=ch_temp, y0=gmin_mean - se, y1=gmin_mean + se, 
           angle=90, code=3, length=0.07, col=growth_T)
  })
  axis(2)
  axis(1, at=seq(17.5, 27.5, by=2.5))
  box()
  legend("topleft", c("Ambient",expression(Ambient~+~3~degree*C)), 
         pch=19, pt.cex=1.1, col=palette(), 
         cex=0.8,
         title="Growth T", bty='n')
  
  
}



figure_hakea_gmin <- function(lopw){
  
  abbrev_hak <- function(x)gsub("hakea","H.",x)
  
  lopw2 <- lopw[order(lopw$gmin.w),]
  lopw2 <- subset(lopw2, !is.na(gmin.d))
  
  bar_cols <- c("lightgrey","black")
  par(mar=c(7,4,1,1), las=3, yaxs="i",
      mgp=c(2.5, 0.5,0), tcl=0.2, cex.axis=0.9, cex.lab=1.1)
  barplot(t(as.matrix(lopw2[,2:3])), beside=T,
          ylim=c(0,1.4),
          col=bar_cols,
          names.arg=abbrev_hak(lopw2$species),
          ylab=expression(g[min]~~(mmol~m^-2~s^-1)))
  box()
  legend("topleft", c("Well-watered","Drought stress"), fill=bar_cols,
         bty='n', cex=0.9)
  
}


