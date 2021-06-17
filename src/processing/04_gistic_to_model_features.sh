#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 2:00:00
#SBATCH --mail-user=<email>
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="gisticfeature"
#SBATCH --output=bbcar/out/gistic_to_model_features.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

python bbcar/src/processing/gistic_to_model_features.py
