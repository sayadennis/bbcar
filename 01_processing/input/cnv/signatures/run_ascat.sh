#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=0
#SBATCH --mem=3G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="ascat%a"
#SBATCH --output=/projects/b1131/saya/bbcar/out/run_ascat_%a.out

module purge all
module load R/4.1.1

#### Create sample names TXT file for job array (RUN THIS IN ADVANCE) ####
# ls /projects/b1131/saya/bbcar/01_alignment/tissue/aligned/*.bam | xargs -n1 basename | tr "_" "\n" | grep "t" > /projects/b1131/saya/bbcar/sample_names_tumor_only.txt

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/bbcar/sample_names_germline.txt

cd /projects/b1131/saya/bbcar/data/02b_cnv/signatures

Rscript --vanilla run_ascat.R
