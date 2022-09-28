#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --gres=gpu:a100:1
#SBATCH -t 24:00:00
#SBATCH -n 12
#SBATCH --mail-user=<email>
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="applythres_bbcar"
#SBATCH --output=bbcar/out/06_freq_based_feature_reduction.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

python bbcar/src/processing/06_freq_based_feature_reduction.py
