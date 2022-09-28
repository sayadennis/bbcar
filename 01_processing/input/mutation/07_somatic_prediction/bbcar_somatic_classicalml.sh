#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="bbcarsom_ml"
#SBATCH --output=bbcar/out/bbcar_somatic_classicalml.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

inputdir='/projects/b1042/lyglab/saya/bbcar/04_ml_features_tmp'
labeldir='/projects/b1042/lyglab/saya/bbcar/04_ml_features_tmp'
ixdir='/projects/b1042/lyglab/saya/bbcar/04_ml_features_tmp/somatic_pred_ix'
outdir='/home/srd6051/bbcar/out'
modeldir='/home/srd6051/bbcar/out'

python classical_ml/ClassicalML/run_classical_ml.py \
    --input $inputdir/input.csv \
    --label $labeldir/target.csv \
    --outfn $outdir/bbcar_somatic_prediction_classicalml.csv \
    --indexdir $ixdir \
    --scoring roc_auc \
    --savemodel $modeldir
#
