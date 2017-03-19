
to.pdf(figure_sim(), "output/figure1.pdf", height=6, width=6)
to.pdf(figure_g0g1_cor(lin2015), "output/figure2.pdf", height=4, width=8)
to.pdf(figure_R2g0(lin2015coef, miner), "output/figure3.pdf", height=4, width=8)

to.pdf(figure_gmin_review_3(kerst2), "output/figure4.pdf", height=4, width=6)
#to.pdf(figure_amings_lin2015), "output/figure5.pdf", height=5, width=5)

#to.pdf(figure_wtc4_gmin(wtc4gmin, wtc4gdark), "output/figure5.pdf", height=4, width=8)
to.pdf(figure_wtc4_gmin_2(wtc4gmin), "output/figure5.pdf", height=4, width=6)


convert_all_pdf("output")

