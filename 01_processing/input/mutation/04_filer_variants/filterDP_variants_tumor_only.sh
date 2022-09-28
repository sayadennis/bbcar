#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --array=0-201
#SBATCH -n 1
#SBATCH -t 4:00:00
#SBATCH --mem=3G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="dpfilter_to_%a"
#SBATCH --output=bbcar/out/DPfilter_tumor_only_%a.out

cd /projects/b1131/saya/bbcar/

## Load GATK 
module purge all
module load bcftools/1.10.1

#### Create sample names TXT file for job array (RUN THIS IN ADVANCE) ####
# ls /projects/b1131/saya/bbcar/02_variant_calls/tumor_only/*.vcf > /projects/b1131/saya/bbcar/filtered_vcf_names_tumor_only.txt

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < sample_names_tumor_only.txt

## Set input and output directories
din='/projects/b1131/saya/bbcar/02_variant_calls/tumor_only'
dout='/projects/b1131/saya/bbcar/02_variant_calls/tumor_only'

##########################
#### From generic PON ####
##########################

bcftools filter --include "INFO/DP>=20" --output-type v --output $dout/${input_args[$SLURM_ARRAY_TASK_ID]}_DPfiltered.vcf $din/${input_args[$SLURM_ARRAY_TASK_ID]}_filtered.vcf

########################
#### From BBCAR PON ####
########################

bcftools filter --include "INFO/DP>=20" --output-type v --output $dout/${input_args[$SLURM_ARRAY_TASK_ID]}_DPfiltered_bbcarpon.vcf $din/${input_args[$SLURM_ARRAY_TASK_ID]}_filtered_bbcarpon.vcf
