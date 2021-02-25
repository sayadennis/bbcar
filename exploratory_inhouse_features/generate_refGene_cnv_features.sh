#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -n 4
#SBATCH -t 12:00:00
#SBATCH --mail-user=<email>
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="refgenecnv"
#SBATCH --output=~/bbcar_project/out/generate_refGene_cnv_features.out

module load python/anaconda3

python bbcar_project/scripts/generate_refGene_cnv_features.py
