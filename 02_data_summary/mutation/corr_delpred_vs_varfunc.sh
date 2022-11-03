#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=5G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="corr"
#SBATCH --output=/projects/b1131/saya/bbcar/out/corr_delpred_vs_varfunc.out

module purge all
module load python-miniconda3/4.12.0
source activate bbcarenv

cd /home/srd6051/bbcar/repo/02_data_summary/mutation

python corr_delpred_vs_varfunc.py
