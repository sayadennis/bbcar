#!/bin/bash
#SBATCH -A p30791
#SBATCH -p normal
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=16G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="cnvsigpro"
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/run_SigProfilerExtractor_cnv_matched.out

module purge all
module load python-miniconda3/4.12.0
source activate SigProfilerExtractor

cd /projects/b1131/saya/new_bbcar/data/02b_cnv/

python ~/bbcar/repo/01_processing/cnv/signatures/05_run_SigProfilerExtractor.py
