#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -n 1
#SBATCH -t 4:00:00
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="cp_raw"
#SBATCH --output=bbcar/out/01_copy_raw_data.sh

cd /projects/b1131/saya/bbcar/

## Copy tissue sequencing raw data 
cp -r /projects/b1122/Zexian/Alignment/BBCAR/RAW_data/* ./raw/tissue/

## Copy germline sequencing raw data 
cp -r /projects/b1122/Zexian/Alignment/Germline_37/RAW_data/* ./raw/germline/
