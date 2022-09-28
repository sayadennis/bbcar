#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --array=0-201
#SBATCH -N 1
#SBATCH -t 48:00:00
#SBATCH --mem=5G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="snpeff_to_\${SLURM_ARRAY_TASK_ID}"
#SBATCH --output=bbcar/out/snpeff_tumor_only_%a.out

cd /projects/b1131/saya/bbcar/

## Load necessary modules 
module purge all
module load java/jdk11.0.10

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/bbcar/sample_names_tumor_only.txt

## Set input and output directories 
din='/projects/b1131/saya/bbcar/02_variant_calls/tumor_only'
dout='/projects/b1131/saya/bbcar/03_annotated_variants/snpeff/tumor_only'
dspf='/projects/b1131/saya/bbcar/tools/snpEff'

java -Xmx5120m -jar ${dspf}/snpEff.jar hg38 ${din}/${input_args[$SLURM_ARRAY_TASK_ID]}_DPfiltered.vcf -stats ${dout}/${input_args[$SLURM_ARRAY_TASK_ID]}_snpEff_genes.txt > ${dout}/${input_args[$SLURM_ARRAY_TASK_ID]}_snpeff.vcf
