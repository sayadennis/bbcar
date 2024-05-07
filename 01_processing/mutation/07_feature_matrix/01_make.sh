#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=5G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="som_feature"
#SBATCH --output=/projects/b1131/saya/bbcar/out/01_make_feature_matrix.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

cd /home/srd6051/bbcar/repo/01_processing/input/mutation/08_feature_matrix/

python 01_make.py
