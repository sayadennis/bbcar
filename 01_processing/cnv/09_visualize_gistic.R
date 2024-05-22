library(maftools)

din <- "/projects/b1131/saya/new_bbcar/data/02b_cnv/09_gistic2_out_conf95"

laml.gistic <- readGistic(gisticDir = dn, isTCGA = FALSE)

fdr <- 0.01

del_cytos <- g[g$Variant_Classification == "Del"]
del_cytos_top20 <- del_cytos[order(del_cytos$qvalues)]$Cytoband[1:20]

amp_cytos <- g[g$Variant_Classification == "Amp"]
amp_cytos_top20 <- amp_cytos[order(amp_cytos$qvalues)]$Cytoband[1:20]

png("gistic_chrom_plot.png", width = 1200, height = 800, res = 150)
gisticChromPlot(
  gistic = laml.gistic,
  fdrCutOff = fdr,
  markBands = union(del_cytos_top20, amp_cytos_top20)
)
dev.off()
