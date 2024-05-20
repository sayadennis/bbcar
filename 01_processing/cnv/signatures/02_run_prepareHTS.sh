#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH --array=0-239
#SBATCH --mem=16G
#SBATCH --job-name=prepHTS
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/run_prepareHTS%a.out

module purge all
module load R/4.1.1

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_all.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

## Obtain IDs of samples that have tissue or germline
IFS=$'\n' read -d '' -r -a tissue < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_tissue.txt
IFS=$'\n' read -d '' -r -a germline < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_germline.txt

## Create output directory
mkdir -p /projects/b1131/saya/new_bbcar/data/02b_cnv/signatures/02_logR_BAF/

if [[ " ${germline[*]} " =~ " ${sampleid} " ]] && [[ " ${tissue[*]} " =~ " ${sampleid} " ]]; then
    Rscript --vanilla ~/bbcar/repo/01_processing/cnv/signatures/02_run_prepareHTS.R $sampleid
fi
