#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=0-201
#SBATCH --mem=12G
#SBATCH --job-name=denoise%a
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/b1131/saya/bbcar/out/denoise_counts%a.out

cd /projects/b1131/saya/bbcar/

## Load GATK 
module purge all 
export PATH="/projects/b1131/saya/bbcar/tools/gatk-4.2.5.0:$PATH"

## Set directories
TISSUE_HDF_DIR="/projects/b1131/saya/bbcar/data/02b_cnv/01_collected_counts/tissue"
PON_DIR="/projects/b1131/saya/bbcar/data/02b_cnv/02_pon"
DENOISED_CTS_DIR="/projects/b1131/saya/bbcar/data/02b_cnv/03_denoised_counts"

mkdir -p $DENOISED_CTS_DIR

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/bbcar/data/sample_ids_all_tissue.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

#### Run for tissue ####
gatk --java-options "-Xmx12g" DenoiseReadCounts \
    -I $TISSUE_HDF_DIR/${sampleid}.counts.hdf5 \
    --count-panel-of-normals ${PON_DIR}/cnvpon.pon.hdf5 \
    --standardized-copy-ratios ${DENOISED_CTS_DIR}/${sampleid}.standardizedCR.tsv \
    --denoised-copy-ratios ${DENOISED_CTS_DIR}/${sampleid}.denoisedCR.tsv;

