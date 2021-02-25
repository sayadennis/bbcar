#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -n 4
#SBATCH -t 2:00:00
#SBATCH --mail-user=<email>
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="targetbbcar"
#SBATCH --output=~/bbcar_project/out/create_target.out

module load python

python bbcar_project/scripts/create_target.py
