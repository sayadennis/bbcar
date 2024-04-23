#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=0-239
#SBATCH --mem=3G
#SBATCH --job-name=plotdenoise%a
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/plot_denoised_copy_ratio%a.out

cd /projects/b1131/saya/new_bbcar/

## Load GATK 
module purge all 
module load R/3.3.1
module load singularity

gatk() {
    singularity exec -B /projects:/projects /projects/p30791/gatk_4.5.0.sif gatk "$@"
}

## Set directories
PON_DIR="/projects/b1131/saya/new_bbcar/data/02b_cnv/02_pon"
DENOISED_CTS_DIR="/projects/b1131/saya/new_bbcar/data/02b_cnv/03_denoised_counts"
PLOT_DIR="/projects/b1131/saya/new_bbcar/plots/cnv/denoised_copy_ratios"
DICT="/projects/b1122/references/GRCH38/v0/Homo_sapiens_assembly38.dict"

mkdir -p $PLOT_DIR

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_tissue.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

gatk PlotDenoisedCopyRatios \
    --standardized-copy-ratios ${DENOISED_CTS_DIR}/${sampleid}.standardizedCR.tsv \
    --denoised-copy-ratios ${DENOISED_CTS_DIR}/${sampleid}.denoisedCR.tsv \
    --sequence-dictionary $DICT \
    --minimum-contig-length 46709983 \
    --output $PLOT_DIR \
    --output-prefix ${sampleid};

