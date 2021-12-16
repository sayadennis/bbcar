#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 00:30:00
#SBATCH -n 4
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="add_baseline"
#SBATCH --output=bbcar/out/bbcar_additional_features_baseline.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

python bbcar/src/modeling/baseline_non_cnv/bbcar_additional_features_baseline.py > bbcar/model_performance/noCNV/baseline_nonCNV_features.txt
