#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=0-239
#SBATCH --mem=96G
#SBATCH --job-name=modelseg%a
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/model_segments%a.out

cd /projects/b1131/saya/new_bbcar/

## Load GATK 
module purge all 
module load singularity

gatk() {
    singularity exec -B /projects:/projects /projects/p30791/gatk_4.5.0.sif gatk "$@"
}

## References etc.
REF="/projects/p30791/hg38_ref/hg38.fa"
DENOISED_CN_DIR="/projects/b1131/saya/new_bbcar/data/02b_cnv/03_denoised_counts"
ALLELIC_CTS_DIR="/projects/b1131/saya/new_bbcar/data/02b_cnv/05_allelic_counts"
CONTIGUOUS_CN_DIR="/projects/b1131/saya/new_bbcar/data/02b_cnv/06_contiguous_cn_segs"

mkdir -p $CONTIGUOUS_CN_DIR/tissue_only/
mkdir -p $CONTIGUOUS_CN_DIR/tissue_normal/

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_tissue.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

## Obtain IDs of samples processed at U Chicago (need to use different interval file) 
IFS=$'\n' read -d '' -r -a uchicago < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_batch_1.txt

## Obtain IDs of samples that have germline
IFS=$'\n' read -d '' -r -a germline < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_germline.txt

#### Run tissue-only for all samples ####
gatk --java-options "-Xmx72g" ModelSegments \
    --denoised-copy-ratios ${DENOISED_CN_DIR}/${sampleid}.denoisedCR.tsv \
    --allelic-counts ${ALLELIC_CTS_DIR}/tissue/${sampleid}.allelicCounts.tsv \
    --output ${CONTIGUOUS_CN_DIR}/tissue_only \
    --output-prefix ${sampleid} \
    --number-of-changepoints-penalty-factor 3.0 \
    --kernel-variance-allele-fraction 0.2 \
    --kernel-variance-copy-ratio 0.2 \
    --kernel-scaling-allele-fraction 0.75 \
    --smoothing-credible-interval-threshold-allele-fraction 4.5 \
    --smoothing-credible-interval-threshold-copy-ratio 4.5

#### Run for tissue-germline pairs if germline is present ####
mkdir -p ${ALLELIC_CTS_DIR}/germline/
if [[ " ${germline[*]} " =~ " ${sampleid} " ]]; then
    gatk --java-options "-Xmx72g" ModelSegments \
        --denoised-copy-ratios ${DENOISED_CN_DIR}/${sampleid}.denoisedCR.tsv \
        --allelic-counts ${ALLELIC_CTS_DIR}/tissue/${sampleid}.allelicCounts.tsv \
        --normal-allelic-counts ${ALLELIC_CTS_DIR}/germline/${sampleid}.allelicCounts.tsv \
        --output ${CONTIGUOUS_CN_DIR}/tissue_normal \
        --output-prefix ${sampleid} \
        --number-of-changepoints-penalty-factor 3.0 \
        --kernel-variance-allele-fraction 0.2 \
        --kernel-variance-copy-ratio 0.2 \
        --kernel-scaling-allele-fraction 0.75 \
        --smoothing-credible-interval-threshold-allele-fraction 4.5 \
        --smoothing-credible-interval-threshold-copy-ratio 4.5
fi

