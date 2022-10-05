#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="bbcarsom_ml"
#SBATCH --output=/projects/b1131/saya/bbcar/out/bbcar_somatic_classicalml.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

cd /projects/b1131/saya/bbcar

inputdir='/projects/b1131/saya/bbcar/data/02a_mutation/04_ml_features'
labeldir='/projects/b1131/saya/bbcar/data/02a_mutation/04_ml_features'
ixdir='/projects/b1131/saya/bbcar/data/02a_mutation/04_ml_features/somatic_pred_ix'
outdir='/projects/b1131/saya/bbcar/model_interpretations'
modeldir='/projects/b1131/saya/bbcar/models'

python ${HOME}/classical-ml/ClassicalML/run_classical_ml.py \
    --input $inputdir/input_matched.csv \
    --label $labeldir/target_matched.csv \
    --outfn $outdir/bbcar_somatic_prediction_classicalml.csv \
    --indexdir $ixdir \
    --scoring roc_auc \
    --savemodel $modeldir
#
