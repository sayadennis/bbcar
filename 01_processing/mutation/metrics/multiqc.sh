#!/bin/bash
#SBATCH -A p30791
#SBATCH -p short
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=12G
#SBATCH --job-name=multiqc
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/multiqc.out


module purge all
module load multiqc/1.2

cd /projects/b1131/saya/new_bbcar/metrics/fastqc_trimmed/

multiqc .

