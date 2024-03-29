#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 2:00:00
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="gistic"
#SBATCH --output=/projects/b1131/saya/bbcar/out/03_run_gistic2_conf90.out

module purge all
module load python-miniconda3/4.12.0
source activate bbcarenv

###############################
#### Process GATK Segments ####
###############################

python 09_process_gatk.py "/projects/b1131/saya/bbcar/data/02b_cnv/07_called_cn_segs/tissue_only"

#####################
#### Run GISTIC2 ####
#####################

din="/projects/b1131/saya/bbcar/data/02b_cnv/07_called_cn_segs/tissue_only"
dout="/projects/b1131/saya/bbcar/data/02b_cnv/09_gistic2_out_conf90"

cd /projects/p30791/GISTIC2/

refgenefile="./refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat"

./gistic2 \
    -b $dout \
    -seg $din/gistic_input_all.tsv \
    -refgene $refgenefile \
    -genegistic 0 \
    -smallmem 0 \
    -brlen 0.98 \
    -conf 0.90 \
    -armpeel 0 \
    -savegene 1 \
    -gcm extreme

###############################
#### Process GATK Segments ####
###############################

python ~/bbcar/repo/01_processing/input/cnv/09_gistic_to_model_features.py

