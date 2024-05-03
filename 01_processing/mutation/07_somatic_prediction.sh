#!/bin/bash
#SBATCH -A p30791
#SBATCH -p long
#SBATCH -t 168:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=48G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=som_ml
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/somatic_prediction.out

module purge all
module load python-miniconda3/4.12.0
source activate classical-ml

cd ${HOME}/classical-ml/ClassicalML/

inputdir='/projects/b1131/saya/new_bbcar/data/02a_mutation/04_ml_features'
labeldir='/projects/b1131/saya/new_bbcar/data/02a_mutation/04_ml_features'
ixdir='/projects/b1131/saya/new_bbcar/data/02a_mutation/04_ml_features/somatic_pred_ix'
outdir='/projects/b1131/saya/new_bbcar/model_interpretations/somatic_prediction'

mkdir -p $outdir

python ${HOME}/classical-ml/ClassicalML/run_classical_ml.py \
    --input $inputdir/input_matched.csv \
    --label $labeldir/target_matched.csv \
    --outdir $outdir \
    --indexdir $ixdir \
    --n_cpu ${SLURM_NTASKS} \
    --scoring roc_auc

cd ${HOME}/bbcar/repo/01_processing/mutation/

python 07_somatic_prediction.py

