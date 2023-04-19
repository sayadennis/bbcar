#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 3:00:00
#SBATCH --mem=12G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=bbcar_freq
#SBATCH --output=/projects/b1131/saya/bbcar/out/03_calculate_bbcar_freq.out

module purge all
module load python-miniconda3/4.12.0
source activate bbcarenv

cd ${HOME}/bbcar/repo/01_processing/input/mutation/exploratory_winham_filter/06_generate_ml_features/

for filter_type in liberal classical strict; do
    python 03_calculate_bbcar_freq.py  $filter_type
done
