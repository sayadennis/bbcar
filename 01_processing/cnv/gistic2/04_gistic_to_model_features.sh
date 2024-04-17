#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 12:00:00
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="gisticfeature"
#SBATCH --output=/projects/b1131/saya/bbcar/out/04_gistic_to_model_features.out

module purge all
module load python-miniconda3/4.12.0
source activate bbcarenv

python ~/bbcar/repo/01_processing/input/cnv/gistic2/04_gistic_to_model_features.py
