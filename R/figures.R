

figure1 <- function(gdfr){
  
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


figure2 <- function(lin2015){
  
  x <- subset(lin2015, fitgroup == "Nicolas Martin-StPaul_Quercus ilex_StPaul_Puechabon")
  fit <- lm(Cond ~ BBopti, data=x)
  
  el <- confidenceEllipse(fit, draw=FALSE)
  set.seed(123)
  b <- bootCase(fit)
  
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




figure3 <- function(lin2015coef, miner){
  
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
  plot_loess(l, add=TRUE, band=FALSE, lwd=2, col="darkgrey")
  
  with(miner,plot(R2, g0, pch=16, ylim=c(-0.1, 0.3), xlim=c(0,1),
                  xlab=expression(R^2),
                  ylab=expression(g[0]~~(mol~m^-2~s^-1))))
  abline(h=0)
  l <- loess(g0 ~ R2, data=miner, span=0.8)
  plot_loess(l, add=TRUE, band=FALSE, lwd=2, col="darkgrey")
  
}


figure4 <- function(lin2015){
  
  par(mar=c(5,5,1,1), cex.axis=0.9, cex.lab=1.1)
  minags <- group_by(lin2015, fitgroup) %>%
    summarize(
      Amin = min(Photo, na.rm=TRUE),
      gsmin = min(Cond, na.rm=TRUE),
      Qrange = max(PARin) - min(PARin),
      Pathway = unique(Pathway)[1]
    ) %>%
    filter(Amin < 5)
  
  with(minags, plot(Amin, gsmin, pch=16,
                    ylim=c(0,0.15),
                    xlab=expression(min~A~~(mu*mol~m^-2~s^-1)),
                    ylab=expression(min~g[s]~(mol~m^-2~s^-1)),
                    panel.first=add_regres_line(lm(gsmin~Amin))))
  
}


figure5 <- function(){
  
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
  
  with(r0, plot(PPFD, ALEAF/GS, type='l', ylim=c(0,80), col=col0, lty=lty0))
  with(r1, lines(PPFD, ALEAF/GS, col=col1, lty=lty1))
  
  with(r0, plot(PPFD, Ci, type='l',
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
                col=col0, lty=lty0,
                ylim=c(0,10)))
  with(r1, lines(Tleaf, ELEAF, col=col1, lty=lty1))
  
  with(r0, plot(VPD, Ci, type='l', 
                col=col0, lty=lty0,
                ylim=c(0,400)))
  with(r1, lines(VPD, Ci, col=col1, lty=lty1))
  
}

