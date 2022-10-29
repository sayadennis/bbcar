#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=3G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="mut_summary"
#SBATCH --output=/projects/b1131/saya/bbcar/out/data_summary_somatic_matrix.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

cd /home/srd6051/bbcar/repo/02_data_summary/mutation

python somatic_matrix.py
