#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 12:00:00
#SBATCH -n 4
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="plotscanmap"
#SBATCH --output=bbcar/out/bbcar_plot_scanmap_training.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

python bbcar/src/modeling/scanmap/bbcar_plot_scanmap_training.py
