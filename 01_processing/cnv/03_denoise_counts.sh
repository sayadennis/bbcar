#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=0-239
#SBATCH --mem=12G
#SBATCH --job-name=denoise%a
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/denoise_counts%a.out

cd /projects/b1131/saya/new_bbcar/

## Load GATK 
module purge all 
module load singularity

gatk() {
    singularity exec -B /projects:/projects /projects/p30791/gatk_4.5.0.sif gatk "$@"
}

## Set directories
TISSUE_HDF_DIR="/projects/b1131/saya/new_bbcar/data/02b_cnv/01_collected_counts/tissue"
PON_DIR="/projects/b1131/saya/new_bbcar/data/02b_cnv/02_pon"
DENOISED_CTS_DIR="/projects/b1131/saya/new_bbcar/data/02b_cnv/03_denoised_counts"

mkdir -p $DENOISED_CTS_DIR

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_tissue.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

#### Run for tissue ####
gatk --java-options "-Xmx12g" DenoiseReadCounts \
    -I $TISSUE_HDF_DIR/${sampleid}.counts.hdf5 \
    --count-panel-of-normals ${PON_DIR}/cnvpon.pon.hdf5 \
    --standardized-copy-ratios ${DENOISED_CTS_DIR}/${sampleid}.standardizedCR.tsv \
    --denoised-copy-ratios ${DENOISED_CTS_DIR}/${sampleid}.denoisedCR.tsv;

