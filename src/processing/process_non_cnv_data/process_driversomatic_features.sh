#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -n 4
#SBATCH --mail-user=<email>
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="driversomatic"
#SBATCH --output=bbcar/out/process_driversomatic_features.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

python bbcar/src/processing/process_non_cnv_data/process_driversomatic_features.py
