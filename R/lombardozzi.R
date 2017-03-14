fn <- "C:\\Users\\30022860\\Google Drive\\Work\\Projects\\WUE Project\\Various\\g0\\Supplementary Table 1 v4_FinalForRevision_MZ.pdf"

#devtools::install_github("ropensci/tabulizerjars", quick=T)
#devtools::install_github("ropensci/tabulizer", quick=T)


library(tabulizer)
out <- extract_tables(fn)

# list, each page is a node

tab <- do.call(c, lapply(out, function(x)x[,3]))

tab <- tab[tab != ""]

gnight <- as.numeric(tab)
gnight <- gnight[!is.na(gnight)]

write.csv(data.frame(gnight=gnight), "data/lombardozzi_gnight.csv", row.names=FALSE)

hist(log10(gnight), axes=FALSE)
axis(2)
library(magicaxis)
magaxis(side=1, unlog=1)

