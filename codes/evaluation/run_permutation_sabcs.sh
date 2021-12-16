#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 6:00:00
#SBATCH -n 4
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="runperm_sabcs"
#SBATCH --output=bbcar/out/run_permutation_sabcs.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

python bbcar/src/evaluation/run_permutation_sabcs.py
