#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="targetbbcar"
#SBATCH --output=/projects/b1131/saya/bbcar/out/05_create_target.out

module purge all
module load python-miniconda3/4.12.0
source activate bbcarenv

python bbcar/src/processing/05_create_target.py
