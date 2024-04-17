#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 3:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="pred_nonmatched"
#SBATCH --output=/projects/b1131/saya/bbcar/out/03_predict_nonmatched.out

module purge all
module load python-miniconda3/4.12.0
source activate classical-ml

cd ${HOME}/bbcar/repo/01_processing/input/mutation/exploratory_winham_filter/07_somatic_prediction/

pon_source='bbcar'

for filter_type in liberal classical strict; do
    dout=/projects/b1131/saya/bbcar/exploratory_winham_filter/${filter_type}/07_predicted_somatic
    mkdir -p $dout
    python 03_predict_nonmatched.py $pon_source $filter_type
done
