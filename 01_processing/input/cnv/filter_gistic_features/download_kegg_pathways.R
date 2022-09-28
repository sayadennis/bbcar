#step 1: install related packages and open them in R:
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("KEGGREST")
# BiocManager::install("EnrichmentBrowser")
library("KEGGREST")
library("EnrichmentBrowser")
library(stringr)
#step2: check and obtain a list of entry identifiers (in this case: sar) and associated definition for a given database or a given set of database entries.
testlist <- keggList("hsa")
#step 3: download the pathways of that organism:
hgpathway <- downloadPathways("hsa")
#step 4: retrieve gene sets for an organism from databases such as GO and KEGG:
hg <- getGenesets(org = "hsa", db = "kegg", cache = TRUE, return.type="list")
# step 4.5 - get a list of names that has word "cancer" in it
cancer_list <- c()
for (item in names(hg)) {
  if (str_detect(item, "[Cc]ancer")) {
    cancer_list <- c(cancer_list, item)
  }
}
#step5: Parse and write the gene sets to a flat text file in GMT format for other pathway enrichment analysis programs (e.g., GSEA):
for (setname in cancer_list) {
  writeGMT(hg["setname"], gmt.file=paste0("/Users/sayadennis/Projects/bbcar_project/kegg_genesets/", setname))
}
