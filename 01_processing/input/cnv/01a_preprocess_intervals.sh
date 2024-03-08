#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=2G
#SBATCH --job-name=procint
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/b1131/saya/bbcar/out/preprocess_intervals.out

## Load GATK 
module purge all 
export PATH="/projects/b1131/saya/bbcar/tools/gatk-4.2.5.0:$PATH"

## Reference genome
REF="/projects/p30791/hg38_ref/hg38.fa"

INT_DIR="/projects/b1131/saya/bbcar/interval_lists"
DICT=/projects/b1122/references/GRCH38/v0/Homo_sapiens_assembly38.dict

## Process V6 interval
gatk BedToIntervalList -I $INT_DIR/SureSelect_v6/S07604514_Regions.bed -O $INT_DIR/SureSelect_v6/hg38.v6.interval_list -SD $DICT

gatk PreprocessIntervals \
    -L "${INT_DIR}/SureSelect_v6/hg38.v6.interval_list" \
    -R $REF \
    --bin-length 0 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O "${INT_DIR}/SureSelect_v6/hg38.v6.preprocessed.interval_list";

## Process V5 interval
gatk BedToIntervalList -I $INT_DIR/SureSelect_v5/S04380110_Regions.bed -O $INT_DIR/SureSelect_v5/hg38.v5.interval_list -SD $DICT

gatk PreprocessIntervals \
    -L "${INT_DIR}/SureSelect_v5/hg38.v5.interval_list" \
    -R $REF \
    --bin-length 0 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O "${INT_DIR}/SureSelect_v5/hg38.v5.preprocessed.interval_list";

