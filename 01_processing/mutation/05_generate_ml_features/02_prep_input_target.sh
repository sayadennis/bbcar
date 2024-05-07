#!/bin/bash
#SBATCH -A p30791
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem=48G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=prepXy
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/prep_input_target.out

module purge all
module load python-miniconda3/4.12.0
source activate bbcarenv

cd ${HOME}/bbcar/repo/01_processing/mutation/06_generate_ml_features/

python 02_prep_input_target.py
