#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=0-201
#SBATCH --mem=72G
#SBATCH --job-name=allelcts%a
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/b1131/saya/bbcar/out/collect_allelic_counts%a.out

cd /projects/b1131/saya/bbcar/

## Load GATK 
module purge all 
export PATH="/projects/b1131/saya/bbcar/tools/gatk-4.2.5.0:$PATH"
module load java/jdk1.8.0_191

## References etc.
REF="/projects/p30791/hg38_ref/hg38.fa"
INT_DIR="/projects/b1131/saya/bbcar/interval_lists"
TISSUE_BAM_DIR="/projects/b1131/saya/bbcar/data/01_alignment/tissue/aligned"
GERMLINE_BAM_DIR="/projects/b1131/saya/bbcar/data/01_alignment/germline/aligned"
ALLELIC_CTS_DIR="/projects/b1131/saya/bbcar/data/02b_cnv/05_allelic_counts"

mkdir -p $ALLELIC_CTS_DIR

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/bbcar/data/sample_ids_all_tissue.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

## Obtain IDs of samples processed at U Chicago (need to use different interval file) 
IFS=$'\n' read -d '' -r -a uchicago < /projects/b1131/saya/bbcar/data/sample_ids_uchicago.txt

## Define interval based on sample ID
if [[ " ${uchicago[*]} " =~ " ${sampleid} " ]]; then
    INTERVAL="${INT_DIR}/SureSelect_v5/hg38.v5.preprocessed.interval_list" # TODO: explore v5 later
else
    INTERVAL="${INT_DIR}/SureSelect_v6/hg38.v6.preprocessed.interval_list"
fi

## Obtain IDs of samples that have germline
IFS=$'\n' read -d '' -r -a germline < /projects/b1131/saya/bbcar/data/sample_ids_all_germline.txt

#### Run for tissue ####
mkdir -p ${ALLELIC_CTS_DIR}/tissue/
gatk --java-options "-Xmx72g" CollectAllelicCounts \
    -L ${INTERVAL} \
    -I ${TISSUE_BAM_DIR}/${sampleid}_bqsr.bam \
    -R ${REF} \
    -O ${ALLELIC_CTS_DIR}/tissue/${sampleid}.allelicCounts.tsv

#### Run for germline if present #### # note: germline is ALWAYS v6 interval!!
mkdir -p ${ALLELIC_CTS_DIR}/germline/
if [[ " ${germline[*]} " =~ " ${sampleid} " ]]; then
gatk --java-options "-Xmx72g" CollectAllelicCounts \
    -L "${INT_DIR}/SureSelect_v6/hg38.v6.preprocessed.interval_list" \
    -I ${GERMLINE_BAM_DIR}/${sampleid}_bqsr.bam \
    -R ${REF} \
    -O ${ALLELIC_CTS_DIR}/germline/${sampleid}.allelicCounts.tsv
fi

