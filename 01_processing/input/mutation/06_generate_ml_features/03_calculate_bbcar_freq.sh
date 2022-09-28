#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem=5G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=bbcar_freq
#SBATCH --output=bbcar/out/03_calculate_bbcar_freq.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

python bbcar/codes/somatic_mutations/02_processing/05_generate_ml_features/03_calculate_bbcar_freq.py
