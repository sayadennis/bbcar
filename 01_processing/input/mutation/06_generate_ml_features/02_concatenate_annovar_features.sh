#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH --mem=5G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=avfts
#SBATCH --output=/projects/b1131/saya/bbcar/out/concatenate_annovar_features.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

cd ${HOME}/bbcar/repo/01_processing/input/mutation/06_generate_ml_features/

#####################
#### Generic PON #### # CURRENTLY THE PYTHON SCRIPT WRITES EVERYTHING TO A SINGLE FILE NAME
#####################

python 02_concatenate_annovar_features.py
