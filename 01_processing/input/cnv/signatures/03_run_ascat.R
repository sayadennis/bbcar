####################################################################################################
## Code taken from instructions in page:
##   https://github.com/VanLoo-lab/ascat/tree/master/ExampleData#processing-targeted-sequencing-data
####################################################################################################

## This needs to be run on an interactive job on HPC because PNG display forwarding is hardcoded \
## and cannot be overwridden by setting arguments. 
## salloc -N 1 -n 1 --account=p30791 --mem=4G --partition=short --time=2:00:00
## module load R/4.1.1
## R
## copy & paste contents of this script in the interactive R session (`Rscript script.R` will not work either)

library(data.table)
library(ASCAT)
source("~/bbcar/repo/01_processing/input/cnv/signatures/custom_ascat_functions.R")

# args <- commandArgs(trailingOnly=TRUE)
# sampleid <- args[1]

germlines <- scan("/projects/b1131/saya/bbcar/data/sample_ids_all_germline.txt", what="", sep="\n") 
tissues <- scan("/projects/b1131/saya/bbcar/data/sample_ids_all_tissue.txt", what="", sep="\n") 

din <- "/projects/b1131/saya/bbcar/data/02b_cnv/signatures/02_logR_BAF"
dout <- "/projects/b1131/saya/bbcar/data/02b_cnv/signatures/03_ASCAT_obj"

sampleids <- scan("/projects/b1131/saya/bbcar/data/sample_ids_all_tissue.txt", what="", sep="\n") 

for (sampleid in sampleids) {
    if ((sampleid %in% germlines) & (sampleid %in% tissues)) {  # case where this is matched tissue-germline 
        ascat.bc = ascat.loadData(
          Tumor_LogR_file = paste0(din, "/tissue_normal/", sampleid, "_tissue_tissueLogR.txt"), 
          Tumor_BAF_file = paste0(din, "/tissue_normal/", sampleid, "_tissue_tissueBAF.txt"), 
          Germline_LogR_file = paste0(din, "/tissue_normal/", sampleid, "_germline_germlineLogR.txt"), 
          Germline_BAF_file = paste0(din, "/tissue_normal/", sampleid, "_germline_germlineBAF.txt"), 
          chrs = paste0('chr',c(1:22, "X")), gender = 'XX', genomeVersion = "hg38", isTargetedSeq=T
        )
    
        # ascat.plotRawData(ascat.bc, img.dir=dout, img.prefix = paste0(sampleid, "_Before_correction_"))
        # ascat.bc = ascat.correctLogR(ascat.bc, GCcontentfile = "GC_file.txt", replictimingfile = "RT_file.txt")
        # ascat.plotRawData(ascat.bc, img.prefix = paste0(sampleid, "_After_correction_"))
        ascat.bc = ascat.aspcf(ascat.bc, penalty=70) # run ASPCF segmentation (?)
        # ascat.plotSegmentedData(ascat.bc, img.dir=dout, img.prefix=paste0(sampleid, "_segment_"))
            
        ascat.output = ascat.runAscat(
          ascat.bc, 
          img.dir=paste0(dout, "/tissue_normal/"), 
          gamma=1, write_segments=TRUE, pdfPlot=TRUE
        )
        
        QC = ascat.metrics(ascat.bc,ascat.output)
        save(ascat.bc, ascat.output, QC, file = paste0(dout, "/", sampleid, "_ASCAT_objects.Rdata"))
    }
}

