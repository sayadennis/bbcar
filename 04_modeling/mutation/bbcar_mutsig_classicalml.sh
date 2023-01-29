#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=6G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="mutsigcml"
#SBATCH --output=/projects/b1131/saya/bbcar/out/bbcar_mutsig_classicalml_all.out

module purge all
module load python-miniconda3/4.12.0
source activate classical-ml

proj_dn="/projects/b1131/saya/bbcar"

inputdir="${proj_dn}/data/02a_mutation/08_feature_matrix/"
labeldir="${proj_dn}/data/clinical/"
outdir="${proj_dn}/model_interpretations/2023_01_mutsig"
ixdir="${proj_dn}/train_test_splits"

#############################
#### Gene-level features ####
#############################

python ${HOME}/classical-ml/ClassicalML/run_classical_ml.py \
    --input $inputdir/signature_per_sample.csv \
    --label $labeldir/bbcar_label_studyid_from_gatk_filenames.csv \
    --outdir $outdir \
    --indexdir $ixdir \
    --scoring roc_auc
#
