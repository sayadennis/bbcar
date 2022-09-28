#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=5G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="eval_sompred"
#SBATCH --output=bbcar/out/evaluate_prediction.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

python bbcar/codes/somatic_mutations/02_processing/06_somatic_prediction/evaluate_prediction.py
