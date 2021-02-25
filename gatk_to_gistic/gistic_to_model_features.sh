#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 2:00:00
#SBATCH --mail-user=<email>
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="gisticfeature"
#SBATCH --output=~/bbcar_project/out/gistic_to_model_features.out

module load python/anaconda3

python bbcar_project/scripts/gistic_to_model_features.py
