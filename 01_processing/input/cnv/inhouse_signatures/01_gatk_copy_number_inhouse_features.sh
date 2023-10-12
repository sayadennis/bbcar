#!/bin/bash
#SBATCH -A p30791
#SBATCH -p short
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=1G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="cnsig"
#SBATCH --output=/projects/b1131/saya/bbcar/out/01_gatk_copy_number_inhouse_features.out

module purge all
module load python-miniconda3/4.12.0
source activate bbcarenv

cd /projects/b1131/saya/bbcar/data/02b_cnv/

python ~/bbcar/repo/01_processing/input/cnv/inhouse_signatures/01_gatk_copy_number_inhouse_features.py
