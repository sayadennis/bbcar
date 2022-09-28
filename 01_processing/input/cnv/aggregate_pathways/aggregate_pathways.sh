#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="aggregate_pathways"
#SBATCH --output=bbcar/out/aggregate_pathways.out

module load R/3.3.1
. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

Rscript --no-save bbcar/src/processing/aggregate_pathways/01_create_hgnc_to_reactome.R
python bbcar/src/processing/aggregate_pathways/02_reactome_subisomorphism.py
python bbcar/src/processing/aggregate_pathways/03_aggregate_cna_to_pathways.py
