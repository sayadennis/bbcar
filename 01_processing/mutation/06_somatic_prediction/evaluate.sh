#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=48G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=someval
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/somatic_prediction_evaluation.out

module purge all
module load python-miniconda3/4.12.0
source activate classical-ml

cd ${HOME}/classical-ml/ClassicalML/

inputdir='/projects/b1131/saya/new_bbcar/data/02a_mutation/04_ml_features'
outputdir='/projects/b1131/saya/new_bbcar/model_interpretations/somatic_prediction_3models'
model_fn='/projects/b1131/saya/new_bbcar/model_interpretations/somatic_prediction_3models/20240523_saved_best_XGB_input_matched.p'

cd ${HOME}/bbcar/repo/01_processing/mutation/06_somatic_prediction/

python evaluate.py \
    $inputdir \
    $outputdir \
    $model_fn

