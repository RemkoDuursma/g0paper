
source("R/load.R")

to.pdf(figure_sim(), "output/figure1.pdf", height=6, width=6)

to.pdf(figure_sureau(planta, plantb), "output/figure2.pdf", height=4, width=9)

to.pdf(figure_gmin_review(gdfr), "output/figure3.pdf", height=5, width=5)

to.pdf(figure_g0g1_cor(lin2015, "Nicolas Martin-StPaul_Quercus ilex_StPaul_Puechabon"), 
       "output/figure4.pdf", height=4, width=8)

to.pdf(figure_R2g0(lin2015coef, miner), "output/figure5.pdf", height=4, width=9)

to.pdf(gmin_3panel(gmindat, cropgmin), "output/figure6.pdf", height=3, width=8)

to.pdf(figure_hakea_gmin(lopw), "output/figure7.pdf", height=5, width=6)

to.pdf(figure_wtc4_gmin_2(wtc4gmin), "output/figure8.pdf", height=4, width=3.5)


#convert_all_pdf("output")


