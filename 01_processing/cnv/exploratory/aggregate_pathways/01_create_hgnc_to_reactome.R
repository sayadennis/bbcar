library(igraph)

## please specify root directory for reactome data
dn.rt <- "/projects/b1122/saya/reactome_pathways"

rerun.rt <- FALSE

if (rerun.rt) {
  library("biomaRt")
  ## Ensembl2Reactome_PE_Pathway.txt is downloaded from Reactome website
  rt <- read.table(
    sprintf("%s/Ensembl2Reactome_PE_Pathway.txt", dn.rt),
    header = FALSE, sep = "\t", quote = "", comment.char = ""
  )
  colnames(rt) <- c(
    "src_db_id", "rt_pe_sid", "rt_pe_nm", "rt_pathway_id",
    "url", "event_nm", "evid_code", "species"
  )

  rt.hs <- rt[rt$species == "Homo sapiens", ]

  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  genes <- rt.hs$src_db_id
  G_list <- unique(getBM(
    filters = "ensembl_gene_id",
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    values = genes, mart = mart
  ))
  rt.hs.gene <- merge(rt.hs, G_list, by.x = "src_db_id", by.y = "ensembl_gene_id")
  save(rt.hs.gene, file = sprintf("%s/Gene2Reactome_PE_Pathway.RData", dn.rt))
} else {
  load(sprintf("%s/Gene2Reactome_PE_Pathway.RData", dn.rt))
}

g2rt <- unique(rt.hs.gene[, c("hgnc_symbol", "rt_pathway_id")])
write.csv(
  g2rt,
  "/projects/b1122/saya/reactome_pathways/gene_to_reactome_pathway.csv",
  row.names = FALSE
)
