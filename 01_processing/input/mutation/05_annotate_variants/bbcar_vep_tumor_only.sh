#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --array=0-201
#SBATCH -N 1
#SBATCH -t 4:00:00
#SBATCH --mem=2G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=vep_to_%a
#SBATCH --output=bbcar/out/vep_tumor_only_%a.out

cd /projects/b1131/saya/bbcar/tools/ensembl-vep

## Load necessary modules 
module purge all
module load perl/5.16

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/bbcar/sample_names_tumor_only.txt

## Set input and output directories 
din='/projects/b1131/saya/bbcar/02_variant_calls/tumor_only'
dout='/projects/b1131/saya/bbcar/03_annotated_variants/vep/tumor_only'
dvep='/projects/b1131/saya/bbcar/tools/ensembl-vep'

${dvep}/vep -i ${din}/${input_args[$SLURM_ARRAY_TASK_ID]}_DPfiltered.vcf -o ${dout}/${input_args[$SLURM_ARRAY_TASK_ID]}_vep.vcf -offline --cache --dir /projects/b1131/saya/bbcar/tools/.vep
