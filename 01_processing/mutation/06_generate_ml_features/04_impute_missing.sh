#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH --mem=16G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=impute_miss
#SBATCH --output=/projects/b1131/saya/bbcar/out/04_impute_missing.out

module purge all
module load python-miniconda3/4.12.0
source activate bbcarenv

cd ${HOME}/bbcar/repo/01_processing/input/mutation/06_generate_ml_features/

python 04_impute_missing.py
