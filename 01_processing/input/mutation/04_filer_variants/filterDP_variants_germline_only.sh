#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --array=0-53
#SBATCH -n 1
#SBATCH -t 4:00:00
#SBATCH --mem=3G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="dpfilter_g%a"
#SBATCH --output=/projects/b1131/saya/bbcar/out/DPfilter_germline_only%a.out

cd /projects/b1131/saya/bbcar/

## Load GATK 
module purge all
module load bcftools/1.10.1

#### Create sample names TXT file for job array (RUN THIS IN ADVANCE) ####
# ls /projects/b1131/saya/bbcar/02_variant_calls/tumor_only/*.vcf > /projects/b1131/saya/bbcar/filtered_vcf_names_tumor_only.txt

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/bbcar/data/02a_mutation/sample_names_germline.txt

## Set input and output directories
din='/projects/b1131/saya/bbcar/data/02a_mutation/02_variant_calls/germline_only'
dout='/projects/b1131/saya/bbcar/data/02a_mutation/02_variant_calls/germline_only'

bcftools filter --include "INFO/DP>=20" --output-type v --output $dout/${input_args[$SLURM_ARRAY_TASK_ID]}_DPfiltered.vcf $din/${input_args[$SLURM_ARRAY_TASK_ID]}_haplotype.vcf
