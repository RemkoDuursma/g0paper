

figure_gmin_review <- function(gdfr){
  
  par(mar=c(4,4,1.5,1), mgp=c(2.5,0.5,0), tcl=-0.2, 
      las=3, yaxs="i",
      cex.lab=1.2)
  plotCIlog(gmin, method, gdfr, 
            ylim=c(-1,2.2),
            label_las=1,
            ylab=expression(Conductance~(mmol~m^-2~s^-1)),
            labels=c(expression(g["cut,isol"]),
                     expression(g["min"]),
                     expression(g["dark"]),
                     expression(g["0"]),
                     expression(g["s,low A"])))
  
}

figure_gmin_review_2 <- function(gdfr){
  
  par(mar=c(4,4,1.5,1), mgp=c(2.5,0.5,0), tcl=-0.2, 
      las=3, yaxs="i",
      cex.lab=1.2)
  plotCI2(gmin, method, gdfr, transform_log10=TRUE,
          label_las=1,
          ylim=c(0,80),
          ylab=expression(Conductance~(mmol~m^-2~s^-1)),
          labels=c(expression(g["cut,isol"]),
                   #expression(g["cut,seal"]),
                   expression(g["min"]),
                   expression(g["dark"]),
                   expression(g["0"]),
                   expression(g["s,low A"])))
  
  
}






figure_gmin_bygroup <- function(gmindat){
  
  par(mar=c(6.5,4,1.5,1), mgp=c(2.5,0.5,0), tcl=-0.2, 
      las=3, yaxs="i",
      cex.lab=1.2, cex.axis=0.8)
  
  plotCI2(gmin, group2, gmindat, transform_log10=TRUE,
          ylab=expression(g[min]~~(mmol~m^-2~s^-1)))
  
}


figure_g0g1_cor <- function(lin2015, group, legend=FALSE){
  
  x <- subset(lin2015, fitgroup == group)
  x$Cond <- 1000*x$Cond
  
  # for bootCase to work
  assign("x",x,envir=.GlobalEnv)
  
  fit <- lm(Cond ~ BBopti, data=x)
  
  el <- confidenceEllipse(fit, draw=FALSE)
  set.seed(123)
  
  b <- bootCase(fit)
  rm(x, envir=.GlobalEnv)
  
  par(mfrow=c(1,2), mar=c(4,4,1,1), mgp=c(2.5,0.5,0), tcl=0.2, las=1,
      cex.lab=1.1, cex.axis=0.9)
  with(x, plot(BBopti, Cond, 
               xlim=c(0,0.04),
               ylim=c(0,150),
               panel.first=add_regres_line(fit),
               xlab=expression(A/(C[a]*sqrt(D))),
               ylab=expression(g[s]~(mmol~m^-2~s^-1)),
               pch=16, cex=0.8, col="dimgrey"))
  box()
  plotlabel("(a)","topleft")
  
  plot(el[,1], el[,2]/1000, type='l', lty=3, 
       xlab=expression(g[0]~~(mmol~m^-2~s^-1)), 
       ylab=expression(g[1]~~(kPa^-0.5)),
       xlim=c(0,25),
       ylim=c(3,5))
  points(b[,1], b[,2]/1000, pch=16, cex=0.3, col="dimgrey")
  points(coef(fit)[1],coef(fit)[2]/1000, pch=19)
  plotlabel("(b)","topleft")
  if(legend)legend("topright", group, bty='n', cex=0.6)
}




figure_R2g0 <- function(lin2015coef, miner){
  
  par(mfrow=c(1,3), mar=c(5,5,1,1), cex.axis=0.9, cex.lab=1.2)
  l <- loess(g0 ~ R2, data=lin2015coef, span=0.8)
  with(lin2015coef, {
    plot(R2, g0, pch=16,
         xlab=expression(R^2),
         ylab=expression(g[0]~~(mmol~m^-2~s^-1)),
         ylim=c(-100, 300),
         xlim=c(0,1),
         panel.first={
           plot_loess(l, add=TRUE, band=TRUE, lwd=2, col="darkgrey")
           panel.last=segments(x0=R2, x1=R2, y0=g0_lci, y1=g0_uci, col="dimgrey")
         }
         
    )
  })
  abline(h=0)
  
  
  legend("topright", "Lin et al. 2015", bty='n', cex=0.8)
  plotlabel("(a)", "topleft")
  
  l <- loess(g0 ~ R2, data=miner, span=0.8)
  with(miner,plot(R2, g0, pch=16, ylim=c(-100, 300), xlim=c(0,1),
                  xlab=expression(R^2),
                  ylab=expression(g[0]~~(mmol~m^-2~s^-1)),
                  panel.first=plot_loess(l, add=TRUE, band=TRUE, lwd=2, col="darkgrey")))
  abline(h=0)
  
  legend("topright", "Miner et al. 2017", bty='n', cex=0.8)
  plotlabel("(b)", "topleft")
  
  l <- loess(g0_se ~ cv_bbopti, data=lin2015coef, span=0.9)
  with(lin2015coef, plot(cv_bbopti, g0_se, pch=16, ylim=c(0,60),
                         xlab=expression("CV of"~A[n]/(C[a] * sqrt(VPD))),
                         ylab=expression("SE of"~g[0]~~(mmol~m^-2~s^-1)),
                         panel.first=plot_loess(l, add=TRUE, band=TRUE, lwd=2, col="darkgrey")
                         ))
  plotlabel("(c)", "topleft")
  
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
  g0_1 <- 0.01
  g0_2 <- 0.03
  
  lty0 <- 1
  lty1 <- 5
  lty2 <- 4
  col0 <- "black"
  col1 <- "black"
  col2 <- "black"
  
  r0 <- Photosyn(PPFD=seq(20,1000, length=101),
                 g1=g1, g0=0, Vcmax=70, Jmax=140,
                 VPD=2, Tleaf=20, Ca=400)
  r1 <- Photosyn(PPFD=seq(20,1000, length=101),
                 g1=g1, g0=g0_1, Vcmax=70, Jmax=140,
                 VPD=2, Tleaf=20, Ca=400)
  r2 <- Photosyn(PPFD=seq(20,1000, length=101),
                 g1=g1, g0=g0_2, Vcmax=70, Jmax=140,
                 VPD=2, Tleaf=20, Ca=400)
  
  with(r0, plot(PPFD, ALEAF/GS, type='l', ylim=c(0,80), col=col0, lty=lty0,
                xlab=expression(PPFD~~(mu*mol~m^-2~s^-1)),
                ylab=expression(A[n]/g[s]~~(mu*mol~mol^-1))
                ))
  with(r1, lines(PPFD, ALEAF/GS, col=col1, lty=lty1))
  with(r2, lines(PPFD, ALEAF/GS, col=col2, lty=lty2))
  plotlabel("(a)", "topright")
  legend("bottomright", c(expression(g[0] == 0), 
                          expression(g[0] == 0.01),
                          expression(g[0] == 0.03)),
         lty=c(lty0, lty1, lty2), bty='n')
  
  with(r0, plot(PPFD, Ci, type='l',
                xlab=expression(PPFD~~(mu*mol~m^-2~s^-1)),
                ylab=expression(C[i]~~(mu*mol~mol^-1)),
                col=col0, lty=lty0,
                ylim=c(250,400)))
  with(r1, lines(PPFD, Ci, col=col1, lty=lty1))
  with(r2, lines(PPFD, Ci, col=col2, lty=lty2))
  plotlabel("(b)", "topright")
  
  tleafs <- seq(10,45,length=101)
  vpds <- 0.000605 * tleafs^2.39
  t0 <- Photosyn(PPFD=1500,
                 g1=g1, g0=0, Vcmax=70, Jmax=140,
                 VPD=vpds, Tleaf=tleafs, Ca=400)
  t1 <- Photosyn(PPFD=1500,
                 g1=g1, g0=g0_1, Vcmax=70, Jmax=140,
                 VPD=vpds, Tleaf=tleafs, Ca=400)
  t2 <- Photosyn(PPFD=1500,
                 g1=g1, g0=g0_2, Vcmax=70, Jmax=140,
                 VPD=vpds, Tleaf=tleafs, Ca=400)
  
  with(t0, plot(Tleaf, 1000*GS, type='l', 
                xlab=expression(T[leaf]~~(degree*C)),
                ylab=expression(g[S]~~(mmol~m^-2~s^-1)),
                col=col0, lty=lty0,
                ylim=c(0,400)))
  with(t1, lines(Tleaf, GS, col=col1, lty=lty1))
  with(t2, lines(Tleaf, GS, col=col2, lty=lty2))
  plotlabel("(c)", "topright")
  
  with(t0, plot(VPD, Ci, type='l', 
                xlab="VPD (kPa)",
                ylab=expression(C[i]~~(mu*mol~mol^-1)),
                col=col0, lty=lty0,
                ylim=c(0,400)))
  with(t1, lines(VPD, Ci, col=col1, lty=lty1))
  with(t2, lines(VPD, Ci, col=col2, lty=lty2))
  plotlabel("(d)", "topright")
  
return(invisible(list(r0=r0, r1=r1, r2=r2, t0=t0, t1=t1, t2=t2)))
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
    arrows(x0=ch_temp, x1=ch_temp, y0=gmin_mean - se, y1=gmin_mean + se, 
           angle=90, code=3, length=0.07, col=growth_T)
  })
  axis(2)
  axis(1, at=seq(17.5, 27.5, by=2.5))
  box()
  legend("topleft", c("Ambient",expression(Ambient~+~3~degree*C)), 
         pch=19, pt.cex=1.1, col=palette(), 
         title="Growth T", bty='n')
  
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
    arrows(x0=xv, x1=xv, y0=gdark_mean - se, y1=gdark_mean + se, 
           angle=90, code=3, length=0.07, col=growth_T)
  })
  axis(2)
  par(mgp=c(2.5,1.5, 0))
  axis(1, at=c(mean(xv[1:2]), mean(xv[3:4])), labels=c("Lower\nsurface", "Upper\nsurface"))
  box()
  
}



figure_wtc4_gmin_2 <- function(wtc4gmin){
  
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
  m <- t(as.matrix(lopw2[,c("gmin.w","gmin.d")]))
  b <- barplot(m, beside=T,
          ylim=c(0,20),
          col=bar_cols,
          names.arg=abbrev_hak(lopw2$species),
          ylab=expression(g[min]~~(mmol~m^-2~s^-1)))
  box()
  legend("topleft", c("Well-watered","Drought stress"), fill=bar_cols,
         bty='n', cex=0.9)
  y <- m[1,]
  x <- colMeans(b)
  ch <- toupper(substr(as.character(lopw2$leaf.form),1,1))
  text(x,y,ch, pos=3, font=3, cex=0.7)
}


gmin_loghist <- function(gmindat){
  v <- log10(gmindat$gmin)
  par(yaxs="i", mgp=c(3,1,0), cex.lab=1.2)
  hist(v, breaks=30, main="", freq=FALSE,
       #ylim=c(0,0.15),
       xlim=c(-1,2),
       col="lightgrey",
       axes=FALSE,
       xlab=expression(g[min]~~(mmol~m^-2~s^-1)))
  magaxis(side=1, unlog=1, usepar=TRUE)
  axis(2)
  curve(dnorm(x, mean=mean(v), sd=sd(v)), add=TRUE,
        n=101)
  
}


gmin_by_order <- function(gmindat){
  
  logmean <- function(x)mean(log(x), na.rm=TRUE)
  
  gmindat3 <- mutate(gmindat, 
                     Order = fct_lump(Order, n=10),
                     Order = reorder(Order, gmin, logmean)) %>%
    filter(Order != "Other") %>% droplevels
  
  plotCI2(gmin, Order, gmindat3, transform_log10=TRUE,
          jit=0.5, datacex=0.5, datacol="grey",
          ylim=c(0,25),
          ylab=expression(g[min]~(mmol~m^-2~s^-1)))

}

gmin_3panel <- function(gmindat, cropgmin){
  par(mfrow=c(1,3), mar=c(6,5,1,1), cex.lab=1.1, mgp=c(2.5, 0.5,0), tcl=-0.2)
  gmin_loghist(gmindat)
  plotlabel("(a)", "topleft")
  par(mar=c(6,5,2,2), las=2, cex.lab=1.1, mgp=c(2.5, 0.5,0), tcl=-0.2)
  gmin_by_order(gmindat)
  plotlabel("(b)", "topleft")
  par(las=2)
  figure_crop_genotypes(cropgmin)
  plotlabel("(c)", "topleft")
}



gmin_by_family <- function(gmindat){
  
  par(mar=c(7,4,2,2))
  gmindat4 <- mutate(gmindat, 
                     Family = fct_lump(Family, n=10),
                     Family = reorder(Family, gmin, median)) %>%
    filter(Family != "Other") %>% droplevels
  
  plotCI2(gmin, Family, gmindat4, transform_log10=TRUE,
          ylab=expression(Conductance~(mmol~m^-2~s^-1)))

}


figure_crop_genotypes <- function(cropgmin){
  
  set.seed(1)
  
  # average by genotypes
  cropgmin2 <- summaryBy(. ~ genotype +crop, data=cropgmin, FUN=mean, na.rm=TRUE, keep.names=TRUE)
  
  d2 <- summaryBy(gmin ~ crop, FUN=c(mean, min, max, length), data=cropgmin2)
  d2 <- d2[order(d2$gmin.mean),]
  n <- nrow(d2)
  
  cropgmin2$crop <- with(cropgmin2, reorder(crop, gmin, mean, na.rm=TRUE))
  cropl <- split(cropgmin2, cropgmin2$crop)
  
  plot(1, type='n', xlim=c(0.4, n + 0.6), ylim=c(0,35), axes=FALSE,
       xlab="",
       ylab=expression(g[min]~~(mmol~m^-2~s^-1)))
  axis(1, at=1:n, labels=d2$crop)
  axis(2)
  box()
  for(i in 1:n){
    with(cropl[[i]], points(jitter(rep(i, nrow(cropl[[i]]))), 
                            gmin, pch=16, col="grey", cex=0.6))
    segments(x0=i,x1=i, y0=d2$gmin.min[i], y1=d2$gmin.max[i])
  }
  points(1:n, d2$gmin.mean, pch=19)
  
  mtext(side=3, at=1:n, line=0.4, text=d2$gmin.length, las=1, cex=0.7)
  mtext(side=3, at=0.25, line=0.4, text="n = ", las=1, cex=0.7)
  
}



figure_sureau <- function(planta, plantb){
  
  #windows(9,4)
  par(mfrow=c(1,3), mar=c(4,4,1,1), mgp=c(2.5, 0.5, 0), tcl=0.2, cex.lab=1.2)
  
  cols <- c(rgb(230, 159, 0, maxColorValue = 255), rgb(0, 114, 178, maxColorValue = 255))
  
  # panel a
  with(planta, plot(days, REW, lwd=2, ylim=c(0,1), xlim=c(0,100), type='l', xlab="Days", col=cols[1]))
  with(plantb, lines(days, REW, lwd=2, ylim=c(0,1), xlim=c(0,100), col=cols[2]))
  plotlabel("(a)", "bottomright")
  legend("topright", c("3.3", "1.6"), title=expression(g[min]), bty='n', lty=1, col=cols, lwd=2, cex=1.2)
  
  # panel b
  with(planta, plot(days, Pmin, ylim=c(-7,0), xlim=c(0,100), 
                    col=cols[1],
                    ylab= "Water potential (MPa)",
                    type='l', lwd=2, lty=1, xlab="Days"))
  with(planta, lines(days, Psoil, lwd=2, col=cols[1], lty=5))
  with(plantb, lines(days, Pmin, type='l', lwd=2, lty=1, col=cols[2]))
  with(plantb, lines(days, Psoil, lwd=2, col=cols[2], lty=5))
  plotlabel("(b)", "bottomright")
  legend("topright", c("leaf","soil"), lty=c(1, 5), lwd=2, bty='n', cex=1.2)
  
  # panel c
  with(planta, 
       plot(days, PLC*100, xlim=c(0,100), 
            ylab="PLC (%)", type='l', lwd=2, 
            xlab="Days", col=cols[1],
            panel.first= abline(h = 88, lty=2, col="darkgrey"))
  )
  with(plantb, lines(days, PLC*100, xlim=c(0,100), type='l', lwd=2, col=cols[2]))
  plotlabel("(c)", "bottomright")
  
}



figure_leafage <- function(gminall){
  
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
  
  ggplot(dat, aes(x=leafage, y=gmin)) +
    geom_line(aes(group=species)) + 
    geom_point() + 
    geom_hline(yintercept=0, alpha=0) + 
    facet_wrap(~datasource, drop=TRUE, scales= "free") +
    theme_bw() +
    scale_colour_manual(values = rainbow(12)) +
    labs(x="", y=expression(g[min]~~(mmol~m^-2~s^-1)))
  
}

