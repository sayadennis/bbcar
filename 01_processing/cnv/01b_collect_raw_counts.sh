#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=0-239
#SBATCH --mem=1G
#SBATCH --job-name=colraw
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/collect_raw_counts%a.out

cd /projects/b1131/saya/new_bbcar/

## Load GATK 
module purge all 
module load singularity

gatk() {
    singularity exec -B /projects:/projects /projects/p30791/gatk_4.5.0.sif gatk "$@"
}

## References etc.
REF="/projects/p30791/hg38_ref/hg38.fa"
INT_DIR="/projects/b1131/saya/bbcar/interval_lists"
TISSUE_BAM_DIR="/projects/b1131/saya/new_bbcar/data/01_alignment/tissue/aligned"
GERMLINE_BAM_DIR="/projects/b1131/saya/new_bbcar/data/01_alignment/germline/aligned"
COLLECTED_CTS_DIR="/projects/b1131/saya/new_bbcar/data/02b_cnv/01_collected_counts"

mkdir -p $COLLECTED_CTS_DIR

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_all.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

## Obtain IDs of samples processed at U Chicago (need to use different interval file) 
IFS=$'\n' read -d '' -r -a uchicago < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_uchicago.txt

## Define interval based on sample ID
if [[ " ${uchicago[*]} " =~ " ${sampleid} " ]]; then
    INTERVAL="${INT_DIR}/SureSelect_v6/hg38.v6.preprocessed.interval_list" # TODO: explore v5 later
else
    INTERVAL="${INT_DIR}/SureSelect_v6/hg38.v6.preprocessed.interval_list"
fi

## Obtain IDs of samples that have germline
IFS=$'\n' read -d '' -r -a tissue < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_tissue.txt
IFS=$'\n' read -d '' -r -a germline < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_germline.txt

#### Run for tissue if present ####
mkdir -p ${COLLECTED_CTS_DIR}/tissue/
if [[ " ${tissue[*]} " =~ " ${sampleid} " ]]; then
    gatk CollectReadCounts \
        -I ${TISSUE_BAM_DIR}/${sampleid}_bqsr.bam \
        -L ${INTERVAL} \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O ${COLLECTED_CTS_DIR}/tissue/${sampleid}.counts.hdf5;
fi

#### Run for germline if present #### # note: germline is ALWAYS v6 interval!!
mkdir -p ${COLLECTED_CTS_DIR}/germline/
if [[ " ${germline[*]} " =~ " ${sampleid} " ]]; then
    gatk CollectReadCounts \
        -I ${GERMLINE_BAM_DIR}/${sampleid}_bqsr.bam \
        -L "${INT_DIR}/SureSelect_v6/hg38.v6.preprocessed.interval_list" \
        --interval-merging-rule OVERLAPPING_ONLY \
        -O ${COLLECTED_CTS_DIR}/germline/${sampleid}.counts.hdf5;
fi

