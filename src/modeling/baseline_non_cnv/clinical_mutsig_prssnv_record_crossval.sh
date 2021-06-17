#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -n 8
#SBATCH --mem=0
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="addfeatcv"
#SBATCH --output=/home/srd6051/bbcar_project/out/210318_bbcar_clinical_mutsig_prssnv_record_crossval.out

module load python/anaconda3

python bbcar_project/scripts/210318_bbcar_clinical_mutsig_prssnv_record_crossval.py
