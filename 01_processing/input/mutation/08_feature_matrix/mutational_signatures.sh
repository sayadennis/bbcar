#!/bin/bash
#SBATCH -A p30791
#SBATCH -p short
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="sub_mutsig"
#SBATCH --output=/projects/b1131/saya/bbcar/out/mutational_signature_subset.out

module purge all
module load python-miniconda3/4.12.0
source activate SigProfilerExtractor

cd /projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix/

python ~/bbcar/repo/01_processing/input/mutation/08_feature_matrix/mutational_signatures.py
