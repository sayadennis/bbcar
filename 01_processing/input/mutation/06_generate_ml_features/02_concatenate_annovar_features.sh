#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH --mem=5G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=avfts
#SBATCH --output=bbcar/out/concatenate_annovar_features.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

#####################
#### Generic PON #### # CURRENTLY THE PYTHON SCRIPT WRITES EVERYTHING TO A SINGLE FILE NAME
#####################

python bbcar/codes/somatic_mutations/02_processing/05_generate_ml_features/concatenate_annovar_features.py
