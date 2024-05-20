library(ASCAT)
source("~/bbcar/repo/01_processing/cnv/signatures/custom_ascat_functions.R")

args <- commandArgs(trailingOnly = TRUE)
sampleid <- args[1]

tissue_name <- paste0(sampleid, "_tissue")
germline_name <- paste0(sampleid, "_germline")

din <- "/projects/b1131/saya/new_bbcar/data/02b_cnv/signatures/01_alleleCounts"
dout <- "/projects/b1131/saya/new_bbcar/data/02b_cnv/signatures/02_logR_BAF"

if (!file.exists(paste0(dout, "/tissue_normal"))) {
  dir.create(paste0(dout, "/tissue_normal"))
}
if (!file.exists(paste0(dout, "/tissue_only"))) {
  dir.create(paste0(dout, "/tissue_only"))
}

# Identify WES interval BED file
uchicago_samples <- scan(
  "/projects/b1131/saya/bbcar/data/sample_ids_uchicago.txt",
  what = "",
  sep = "\n"
)

int_dir <- "/projects/b1122/gannon/bbcar/RAW_data/int_lst"

if (sampleid %in% uchicago_samples) {
  int_file <- paste0(int_dir, "/SureSelect_v5/hg38/hg38.preprocessed.interval_list")
} else {
  int_file <- paste0(int_dir, "/SureSelect_v6/hg38.preprocessed.interval_list")
}

# Obtain BAF and LogR from the raw allele counts
custom_getBAFsAndLogRs(
  tissue_name = tissue_name,
  germline_name = germline_name,
  dout = dout,
  tissueAlleleCountsFile.prefix = paste0(din, "/", tissue_name, "_alleleFrequencies"),
  normalAlleleCountsFile.prefix = paste0(din, "/", germline_name, "_alleleFrequencies"),
  alleles.prefix = "/projects/b1131/saya/bbcar/tools/Battenberg/battenberg_1000genomesloci2012_v3/1000genomesAlleles2012_chr", # nolint
  gender = "XX",
  genomeVersion = "hg38",
  chrom_names = c(1:22),
  minCounts = 10,
  BED_file = int_file,
  probloci_file = "/projects/b1131/saya/bbcar/tools/Battenberg/probloci_270415.txt",
  seed = as.integer(Sys.time())
) # nolint
# traceback()
options(error = traceback)

# # Synchronise all information
# ascat.synchroniseFiles(samplename=tissue_name,
#                         tumourLogR_file=tissueLogR_file,
#                         tumourBAF_file=tissueBAF_file,
#                         normalLogR_file=germlineLogR_file,
#                         normalBAF_file=germlineBAF_file
# )
