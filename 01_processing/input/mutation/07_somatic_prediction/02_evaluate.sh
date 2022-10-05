#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=5G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="eval_sompred"
#SBATCH --output=/projects/b1131/saya/bbcar/out/02_evaluate_somatic_pred.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

cd ${HOME}/bbcar/repo/01_processing/input/mutation/07_somatic_prediction/

python 02_evaluate.py
