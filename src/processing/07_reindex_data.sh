#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-gpu
#SBATCH -p genomics
#SBATCH -t 24:00:00
#SBATCH -n 12
#SBATCH --mail-user=<email>
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="reindex_bbcar"
#SBATCH --output=bbcar/out/07_reindex_data.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

python bbcar/src/processing/07_reindex_data.py
