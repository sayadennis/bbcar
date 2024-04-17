#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 4:00:00
#SBATCH --mem=8G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=input_target
#SBATCH --output=/projects/b1131/saya/bbcar/out/05_prep_input_target.out

module purge all
module load python-miniconda3/4.12.0
source activate bbcarenv

cd ${HOME}/bbcar/repo/01_processing/input/mutation/06_generate_ml_features/

python 05_prep_input_target.py
