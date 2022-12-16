####################################################################################################
## Code taken from instructions in page:
##   https://github.com/VanLoo-lab/ascat/tree/master/ExampleData#processing-targeted-sequencing-data
####################################################################################################

library(ASCAT)

ASCATobj <- ascat.loadData(
  Tumor_LogR_file, 
  Tumor_BAF_file, 
  Germline_LogR_file = NULL, 
  Germline_BAF_file = NULL, 
  chrs = c(1:22,"X","Y"), 
  gender = NULL, 
  sexchromosomes = c("X","Y"), 
  genomeVersion=NULL, 
  isTargetedSeq=F
)

germ_geno <- ascat.predictGermlineGenotypes(
  ASCATobj, 
  platform = "AffySNP6", 
  img.dir=".", 
  img.prefix=""
)

saveRDS(germ_geno, "/projects/b1131/saya/bbcar/data/ascat/test_germ_geno.RDS")

ascat.prepareTargetedSeq(
  Worksheet = "myWorksheet.tsv", # A tab-separated file with specific information. Check format using ?ascat.prepareTargetedSeq
  alleles.prefix = "G1000_alleles_hg19_chr",
  BED_file = "my_targeted_design.bed",
  allelecounter_exe = "/PATH/TO/allelecounter",
  genomeVersion = "hg19",
  nthreads = 8)

ascat.prepareHTS(
  tumourseqfile = "Tumour.bam",
  normalseqfile = "Normal.bam",
  tumourname = "Tumour_name",
  normalname = "Normal_name",
  allelecounter_exe = "/PATH/TO/allelecounter",
  alleles.prefix = "./alleleData/Cleaned/alleleData_chr",
  loci.prefix = "./alleleData/Cleaned/loci_chr",
  gender = "XX",
  genomeVersion = "hg19",
  nthreads = 8,
  tumourLogR_file = "Tumor_LogR.txt",
  tumourBAF_file = "Tumor_BAF.txt",
  normalLogR_file = "Germline_LogR.txt",
  normalBAF_file = "Germline_BAF.txt")
  
ascat.bc = ascat.loadData(Tumor_LogR_file = "Tumor_LogR.txt", Tumor_BAF_file = "Tumor_BAF.txt", Germline_LogR_file = "Germline_LogR.txt", Germline_BAF_file = "Germline_BAF.txt", gender = 'XX', genomeVersion = "hg19", isTargetedSeq=T)
ascat.plotRawData(ascat.bc, img.prefix = "Before_correction_")
ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "GC_file.txt", replictimingfile = "RT_file.txt")
ascat.plotRawData(ascat.bc, img.prefix = "After_correction_")
ascat.bc = ascat.aspcf(ascat.bc, penalty=25)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc, gamma=1, write_segments = T)
QC = ascat.metrics(ascat.bc,ascat.output)
save(ascat.bc, ascat.output, QC, file = 'ASCAT_objects.Rdata')
