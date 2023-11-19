#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH -n 1
#SBATCH --mem=12G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="explore_suphNMF"
#SBATCH --output=/projects/b1131/saya/bbcar/out/supervised_hNMF_mutation_cnv.out

module purge all
module load python-miniconda3/4.12.0
source activate bbcarenv

cd ${HOME}/bbcar/repo/04_modeling/

python supervised_hNMF_mutation_cnv.py

python unsupervised_hNMF_mutation_cnv.py
