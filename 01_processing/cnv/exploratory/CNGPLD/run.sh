#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=4G
#SBATCH --job-name=cngpld
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/b1131/saya/bbcar/out/test_CNGPLD.out

module purge all
module load R/4.1.1

cd /projects/p30791/cngpld/

Rscript --vanilla ~/bbcar/repo/01_processing/input/cnv/exploratory/CNGPLD/run.R

