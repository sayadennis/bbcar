#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=10G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="som_agg"
#SBATCH --output=/projects/b1131/saya/bbcar/out/01_aggregate_feature_matrix.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

cd /home/srd6051/bbcar/repo/01_processing/input/mutation/08_feature_matrix/

python 02_aggregate.py
