#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --array=0-201
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH --mem=1G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="affilter_to_%a"
#SBATCH --output=/projects/b1131/saya/bbcar/out/AFfilter_tumor_only_%a.out

cd /projects/b1131/saya/bbcar/

## Load GATK 
module purge all
module load bcftools/1.10.1

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/bbcar/data/02a_mutation/sample_names_tumor_only.txt

## Set input and output directories
din='/projects/b1131/saya/bbcar/data/02a_mutation/02_variant_calls/tumor_only'
dout='/projects/b1131/saya/bbcar/exploratory_winham_filter'

mkdir -p $dout

#################
#### Liberal ####
#################

mkdir -p ${dout}/liberal/02_variant_calls/tumor_only
bcftools filter --include "FORMAT/AF>=0.05" --output-type v --output ${dout}/liberal/02_variant_calls/tumor_only/${input_args[$SLURM_ARRAY_TASK_ID]}_DPfiltered_liberalAF_bbcarpon.vcf $din/${input_args[$SLURM_ARRAY_TASK_ID]}_DPfiltered_bbcarpon.vcf

###################
#### Classical ####
###################

mkdir -p ${dout}/classical/02_variant_calls/tumor_only
bcftools filter --exclude '(FORMAT/AF<0.05) | (FORMAT/AF<0.10 && REF="C" && ALT="T")' --output-type v --output ${dout}/classical/02_variant_calls/tumor_only/${input_args[$SLURM_ARRAY_TASK_ID]}_DPfiltered_classicalAF_bbcarpon.vcf $din/${input_args[$SLURM_ARRAY_TASK_ID]}_DPfiltered_bbcarpon.vcf

################
#### Strict ####
################

mkdir -p ${dout}/strict/02_variant_calls/tumor_only
bcftools filter --include "FORMAT/AF>=0.10" --output-type v --output ${dout}/strict/02_variant_calls/tumor_only/${input_args[$SLURM_ARRAY_TASK_ID]}_DPfiltered_strictAF_bbcarpon.vcf $din/${input_args[$SLURM_ARRAY_TASK_ID]}_DPfiltered_bbcarpon.vcf
