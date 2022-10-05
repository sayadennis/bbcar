#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=5G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="pred_nonmatched"
#SBATCH --output=bbcar/out/03_predict_nonmatched.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

cd /home/srd6051/bbcar/repo/01_processing/input/mutation/07_somatic_prediction/

python 03_predict_nonmatched.py
