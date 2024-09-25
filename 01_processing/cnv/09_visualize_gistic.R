# R/4.3.0

library(tidyverse)
library(maftools)

fdr <- 0.01
g_thres <- 0.75

is_overlapping <- function(chrom1, start1, end1, chrom2, start2, end2) {
  if (chrom1 == chrom2) {
    ol <- (start1 <= end2) && (start2 <= end1)
  } else {
    ol <- FALSE
  }
  ol
}

cytoband_by_g <- function(g, gis.score, g_thres) {
  ols <- list()

  gis.score <- gis.scores[gis.scores$G_Score > g_thres, ]

  for (i in seq_len(nrow(g))) {
    chrom1 <- g$Chromosome[i]
    start1 <- g$Start_Position[i]
    end1 <- g$End_Position[i]
    for (j in seq_len(nrow(gis.scores))) {
      chrom2 <- as.numeric(gis.scores$Chromosome[j])
      start2 <- as.numeric(gis.scores$Start_Position[j])
      end2 <- as.numeric(gis.scores$End_Position[j])
      ol <- is_overlapping(chrom1, start1, end1, chrom2, start2, end2)
      if (ol) {
        gis.scores[j, "Cytoband"] <- g$Cytoband[i]
        ols <- append(ols, list(c(i, j)))
      }
    }
  }

  cytobands <- unique(gis.scores[gis.scores$G_Score > g_thres, ]$Cytoband)
  cytobands <- cytobands[!sapply(cytobands, is.na)]

  cytobands
}

for (group in c("case", "control")) {
  din <- paste0("/projects/b1131/saya/new_bbcar/data/02b_cnv/09_gistic2_out_conf95_", group)

  laml.gistic <- readGistic(gisticDir = din, isTCGA = FALSE)

  g <- getCytobandSummary(laml.gistic)
  g <- g[qvalues < fdr]

  print(group)
  num_amps <- sum(g$Variant_Classification == "Amp")
  num_dels <- sum(g$Variant_Classification == "Del")
  print(paste0("Amplification: ", num_amps, " // Deletions: ", num_dels))

  # Add genomic coordinates as separate columns
  g <- g %>%
    mutate(
      Chromosome = as.integer(str_replace(str_extract(Wide_Peak_Limits, "^chr\\d+"), "chr", "")),
      Start_Position = as.integer(str_extract(str_extract(Wide_Peak_Limits, ":\\d+-"), "\\d+")),
      End_Position = as.integer(str_extract(str_extract(Wide_Peak_Limits, "-\\d+$"), "\\d+"))
    )

  # ## If I want to color top 20 cytobands ##
  # del_cytos <- g[g$Variant_Classification == "Del"]
  # del_cytos_top20 <- del_cytos[order(del_cytos$qvalues)]$Cytoband[1:20]
  #
  # amp_cytos <- g[g$Variant_Classification == "Amp"]
  # amp_cytos_top20 <- amp_cytos[order(amp_cytos$qvalues)]$Cytoband[1:20]
  # cytobands <- union(del_cytos_top20, amp_cytos_top20)

  ## If I want to mark cytobands by G-score
  gis.scores <- maftools:::transformSegments(segmentedData = laml.gistic@gis.scores, build = "hg38")
  gis.scores$amp <- ifelse(
    test = gis.scores$Variant_Classification == "Del",
    yes = -gis.scores$G_Score,
    no = gis.scores$G_Score
  )
  gis.scores$ystart <- ifelse(
    test = gis.scores$Variant_Classification == "Del",
    yes = -0.01,
    no = 0.01
  )
  gis.scores$Variant_Classification <- ifelse(
    test = as.numeric(gis.scores$fdr) > fdr,
    yes = gis.scores$Variant_Classification,
    no = "neutral"
  )
  gis.scores$Variant_Classification <- factor(
    gis.scores$Variant_Classification,
    levels = c("neutral", "Amp", "Del")
  )

  cytobands <- cytoband_by_g(g, gis.scores, g_thres)

  png(paste0("gistic_chrom_plot_", group, ".png"), width = 1500, height = 500, res = 150)
  gisticChromPlot(
    gistic = laml.gistic,
    fdrCutOff = fdr,
    y_lims = c(-2.0, 2.0),
    markBands = cytobands
  )
  dev.off()
  print(paste0("Finished for ", group, "!"))
}
