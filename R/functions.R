to.pdf <- function(expr, filename, ..., verbose=TRUE) {
  if(!file.exists(dirname(filename)))
    dir.create(dirname(filename), recursive=TRUE)
  if ( verbose )
    cat(sprintf("Creating %s\n", filename))
  pdf(filename, ...)
  on.exit(dev.off())
  eval.parent(substitute(expr))
}

plotfig <- function(...)grid.raster(readPNG(...))
pval <- function(...)formatPval(..., includeEquality=TRUE)


fits_lin2015 <- function(lin2015a){
  
  fits <- lmList(Cond ~ BBopti | fitgroup, data=lin2015a, na.action=na.omit)
  p <- as.data.frame(coef(fits)) %>% 
    rownames_to_column(var = "fitgroup") %>%  
    rename(g0 = `(Intercept)`, g1 = BBopti) %>% 
    mutate(g0 = 1000 * g0,
           R2 = sapply(fits, function(x)summary(x)$adj.r.squared),
           rmse = sapply(fits, function(x)summary(x)$sigma),
           cor = sapply(fits, function(x)summary(x, correlation=T)$correlation[2,1]),
           g0_lci = 1000 * sapply(fits, function(x)confint(x)[1,1]),
           g0_uci = 1000 * sapply(fits, function(x)confint(x)[1,2]),
           g0_se = 1000 * sapply(fits, function(x)summary(x)$coefficients[1,2]),
           g1_se = sapply(fits, function(x)summary(x)$coefficients[2,2])
           ) %>% 
    left_join(dplyr::summarize(lin2015a, 
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


# Simple function for placing labels on a figure.
plotlabel <- function(txt, where, inset=0.08, inset.x=inset, inset.y=inset,...){
  u <- par()$usr
  if(grepl("left",where))x <- u[1] + inset.x*(u[2]-u[1])
  if(grepl("right",where))x <- u[2] - inset.x*(u[2]-u[1])
  if(grepl("bottom",where))y <- u[3] + inset.y*(u[4]-u[3])
  if(grepl("top",where))y <- u[4] - inset.y*(u[4]-u[3])
  
  text(x,y,txt,font=2,...)
}


convert_all_pdf <- function(path, ...){
  o <- getwd()
  on.exit(setwd(o))
  setwd(path)
  pdfs <- dir(pattern="[.]pdf")
  for(i in seq_along(pdfs)){
    convert_pdf_png(pdfs[i])
    message(i)
  }
}

convert_pdf_png <- function(filename, fnout=NULL, overwrite=TRUE, res=600, resize=100){
  if(is.null(fnout))fnout <- gsub("[.]pdf",".png",filename)
  if(!file.exists(fnout) | overwrite){
    cmd <- sprintf("convert -density %s -resize %s%% -colorspace CMYK %s %s",res, resize, filename,fnout) 
    shell(cmd)
  }
}

# Simple function for placing labels on a figure.
plotlabel <- function(txt, where, inset=0.08, inset.x=inset, inset.y=inset,...){
  u <- par()$usr
  if(grepl("left",where))x <- u[1] + inset.x*(u[2]-u[1])
  if(grepl("right",where))x <- u[2] - inset.x*(u[2]-u[1])
  if(grepl("bottom",where))y <- u[3] + inset.y*(u[4]-u[3])
  if(grepl("top",where))y <- u[4] - inset.y*(u[4]-u[3])
  
  text(x,y,txt,font=2,...)
}



plotCI2 <- function(yvar, group, data, transform_log10=FALSE, 
                    ylab=NULL, ylim=NULL, labels=NULL, label_las=3, add_data=TRUE, 
                    jit=0.6, datacex=0.4, datacol="lightgrey", 
                    ...){
  
  if(is.null(ylab))ylab <- substitute(yvar)
  
  data <- as.data.frame(data)
  
  data$Y <- eval(substitute(yvar), data)
  data$G <- as.factor(eval(substitute(group), data))
  
  datl <- split(data, data$G)
  
  if(is.null(labels)){
    labels <- levels(data$G)
  }
  
  if(transform_log10){
    fit <- lm(log10(Y) ~ G - 1, data=data)
    ci <- 10^confint(fit)
    mn <- 10^coef(fit)
  } else {
    fit <- lm(Y ~ G - 1, data=data)
    ci <- confint(fit)
    mn <- coef(fit)
  }
  
  if(is.null(ylim))ylim <- c(0, max(ci)+0.1*max(ci))
  g <- glht(fit, linfct=mcp(G = "Tukey"))
  lets <- cld(g)$mcletters$Letters
  n <- length(mn)
  
  plot(1:n, mn, pch=19, cex=1.1, ylim=ylim, axes=FALSE,
       xlab="",
       panel.first={
         if(add_data)with(data, points(jitter(as.numeric(G), factor=jit), data$Y, pch=16, col=datacol, cex=datacex))
       },
       xlim=c(0.7, n+0.3), ylab=ylab, ...)
  axis(1, at=1:n, labels=labels, las=label_las)
  
  arrows(x0=1:n, x1=1:n, y0=ci[,1], y1=ci[,2], angle=90, length=0.05, code=3)
  text(x=1:n, y=ci[,2], lets, pos=3, cex=0.8)
  par(las=1)
  axis(2)
  box()
  
  tb <- unname(table(data$G))
  mtext(side=3, at=1:n, line=0.4, text=tb, las=1, cex=0.7)
  mtext(side=3, at=0.25, line=0.4, text="n = ", las=1, cex=0.7)
  
}


plotCIlog <- function(yvar, group, data, 
                    ylab=NULL, ylim=NULL, labels=NULL, label_las=3, add_data=TRUE, 
                    jit=0.6, datacex=0.4, datacol="lightgrey", 
                    ...){
  
  if(is.null(ylab))ylab <- substitute(yvar)
  
  data <- as.data.frame(data)
  
  data$Y <- eval(substitute(yvar), data)
  data$G <- as.factor(eval(substitute(group), data))
  
  datl <- split(data, data$G)
  
  if(is.null(labels)){
    labels <- levels(data$G)
  }
  
  fit <- lm(log10(Y) ~ G - 1, data=data)
  ci <- confint(fit)
  mn <- coef(fit)

  if(is.null(ylim))ylim <- c(0, max(ci)+0.1*max(ci))
  g <- glht(fit, linfct=mcp(G = "Tukey"))
  lets <- cld(g)$mcletters$Letters
  
  n <- length(mn)
  
  plot(1:n, mn, pch=19, cex=1.1, ylim=ylim, axes=FALSE,
       xlab="",
       panel.first={
         if(add_data)with(data, points(jitter(as.numeric(G), factor=jit), log10(gmin), pch=16, col=datacol, cex=datacex))
       },
       xlim=c(0.7, n+0.3), ylab=ylab, ...)
  axis(1, at=1:n, labels=labels, las=label_las)
  
  arrows(x0=1:n, x1=1:n, y0=ci[,1], y1=ci[,2], angle=90, length=0.05, code=3)
  text(x=1:n, y=ci[,2], lets, pos=3, cex=0.8)
  par(las=1)
  magaxis(side=2, unlog=2)
  box()
  
  tb <- unname(table(data$G))
  mtext(side=3, at=1:n, line=0.4, text=tb, las=1, cex=0.7)
  mtext(side=3, at=0.25, line=0.4, text="n = ", las=1, cex=0.7)
  
}


as_dataframe_cls <- function(cls){
  
  fun <- function(x){
    if(all(is.na(x))){
      return(NULL)
    } else {
      vec <- c(species = x$name[x$rank == 'species'],
               Family = x$name[x$rank == 'family'],
               Order = x$name[x$rank == 'order'])
      if(!"species" %in% names(vec)){
        return(NULL)
      } else {
        return(vec)
      }
    }
  }
  
  m <- do.call(rbind, lapply(cls, fun))
return(as.data.frame(m, stringsAsFactors = FALSE))
}

