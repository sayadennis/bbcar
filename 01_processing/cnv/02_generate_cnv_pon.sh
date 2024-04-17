#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=24G
#SBATCH --job-name=cnvpon
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/b1131/saya/bbcar/out/generate_cnv_pon.out

cd /projects/b1131/saya/bbcar/

## Load GATK 
module purge all 
export PATH="/projects/b1131/saya/bbcar/tools/gatk-4.2.5.0:$PATH"

## Directories
GERMLINE_HDF_DIR="/projects/b1131/saya/bbcar/data/02b_cnv/01_collected_counts/germline"
PON_DIR="/projects/b1131/saya/bbcar/data/02b_cnv/02_pon"

mkdir -p $PON_DIR

gatk --java-options "-Xmx24g" CreateReadCountPanelOfNormals \
    --minimum-interval-median-percentile 5.0 \
    -O $PON_DIR/cnvpon.pon.hdf5 \
    $(for x in $(ls -1 $GERMLINE_HDF_DIR/*.counts.hdf5); do echo -n "-I $x "; done);

