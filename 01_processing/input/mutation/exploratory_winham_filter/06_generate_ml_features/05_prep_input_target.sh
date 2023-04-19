#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH --mem=6G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=input_target
#SBATCH --output=/projects/b1131/saya/bbcar/out/05_prep_input_target.out

module purge all
module load python-miniconda3/4.12.0
source activate bbcarenv

cd ${HOME}/bbcar/repo/01_processing/input/mutation/exploratory_winham_filter/06_generate_ml_features/

for filter_type in liberal classical strict; do
    mkdir -p /projects/b1131/saya/bbcar/exploratory_winham_filter/${filter_type}/04_ml_features/somatic_pred_ix
    python 05_prep_input_target.py $filter_type
done
