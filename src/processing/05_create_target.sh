#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -n 4
#SBATCH -t 2:00:00
#SBATCH --mail-user=<email>
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="targetbbcar"
#SBATCH --output=bbcar/out/create_target.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

python bbcar/src/processing/create_target.py
