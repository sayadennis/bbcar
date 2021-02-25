#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 24:00:00
#SBATCH --mail-user=<email>
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="gistic"
#SBATCH --output=~/bbcar_project/out/run_gistic2.out

basedir="/projects/b1122/saya/bbcar_project/gistic2_out"
segfile="/projects/b1122/saya/bbcar_project/gistic2_input/combined_gistic_input.tsv"
refgenefile="~/GISTIC2/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat"
# intervals="/projects/b1122/gannon/bbcar/RAW_data/int_lst/SureSelect_v6/hg38.preprocessed.interval_list"

cd /home/srd6051/GISTIC2

./gistic2 -b $basedir -seg $segfile -refgene $refgenefile -genegistic 0 -smallmem 0 -brlen 0.98 -conf 0.75 -armpeel 0 -savegene 1 -gcm extreme
# -mk $intervals 
