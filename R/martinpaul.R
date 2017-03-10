
source("R/load.R")

library(car)

x <- subset(lin2015, fitgroup == "Nicolas Martin-StPaul_Quercus ilex_StPaul_Puechabon")
fit <- lm(Cond ~ BBopti, data=x)

el <- confidenceEllipse(fit, draw=FALSE)
b <- bootCase(fit)


windows(8,4)
par(mfrow=c(1,2), mar=c(4,4,1,1), mgp=c(2.5,0.5,0), tcl=0.2, las=1,
    cex.lab=1.1, cex.axis=0.9)
with(x, plot(BBopti, Cond, xlim=c(0, max(BBopti)), ylim=c(0, 0.15),
             panel.first=add_regres_line(fit),
             xlab=expression(A/(C[a]*sqrt(D))),
             ylab=expression(g[s]~(mol~m^-2~s^-1)),
             pch=16, cex=0.8, col="dimgrey"))
box()

plot(el, type='l', lty=3, 
     xlab=expression(g[0]), ylab=expression(g[1]),
     ylim=c(3,5))
points(b, pch=16, cex=0.3, col="dimgrey")
points(coef(fit)[1],coef(fit)[2], pch=19)




confint(lm(Cond ~ BBopti, data=x))





