#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH --array=0-201
#SBATCH --mem=16G
#SBATCH --job-name=prepHTS%a
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/b1131/saya/bbcar/out/run_prepareHTS%a.out

module purge all
module load R/4.1.1

IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/bbcar/data/sample_ids_all_tissue.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

Rscript --vanilla ~/bbcar/repo/01_processing/input/cnv/signatures/02_run_prepareHTS.R $sampleid
