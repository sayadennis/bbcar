#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 24:00:00
#SBATCH --mail-user=<email>
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="inputgistic"
#SBATCH --output=~/bbcar_project/out/gatk_segment_to_gistic_input.out

module load python/anaconda3

python bbcar_project/scripts/gatk_segment_to_gistic_input.py
