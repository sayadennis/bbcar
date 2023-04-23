#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="eval_sompred"
#SBATCH --output=/projects/b1131/saya/bbcar/out/02_evaluate_somatic_pred.out

module purge all
module load python-miniconda3/4.12.0
source activate classical-ml

cd ${HOME}/bbcar/repo/01_processing/input/mutation/07_somatic_prediction/

# for feature_type in bbcar 1000g; do
#     python 02_evaluate.py $feature_type
# done

python 02_evaluate.py bbcar
