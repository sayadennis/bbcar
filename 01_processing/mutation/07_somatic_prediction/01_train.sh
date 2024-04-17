#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicslong
#SBATCH -t 168:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=18G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="bbcarsom_ml"
#SBATCH --output=/projects/b1131/saya/bbcar/out/bbcar_somatic_classicalml.out

module purge all
module load python-miniconda3/4.12.0
source activate classical-ml

cd /projects/b1131/saya/bbcar

inputdir='/projects/b1131/saya/bbcar/data/02a_mutation/04_ml_features'
labeldir='/projects/b1131/saya/bbcar/data/02a_mutation/04_ml_features'

for feature_type in bbcar 1000g; do
    ixdir=/projects/b1131/saya/bbcar/data/02a_mutation/04_ml_features/somatic_pred_ix/${feature_type}
    outdir=/projects/b1131/saya/bbcar/model_interpretations/${feature_type}
    modeldir=/projects/b1131/saya/bbcar/models/${feature_type}

    mkdir -p $outdir
    mkdir -p $modeldir

    python ${HOME}/classical-ml/ClassicalML/run_classical_ml.py \
        --input $inputdir/input_matched_${feature_type}.csv \
        --label $labeldir/target_matched_${feature_type}.csv \
        --outdir $outdir \
        --indexdir $ixdir \
        --scoring roc_auc
done
