#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="mut_cml"
#SBATCH --output=/projects/b1131/saya/bbcar/out/bbcar_mut_classicalml_all.out

module purge all
module load python-miniconda3/4.12.0
source activate classical-ml

proj_dn="/projects/b1131/saya/bbcar"

inputdir="${proj_dn}/data/02a_mutation/08_feature_matrix/"
labeldir="${proj_dn}/data/clinical/"
outdir="${proj_dn}/model_interpretations"
ixdir="${proj_dn}/train_test_splits"

#############################
#### Gene-level features ####
#############################

python ${HOME}/classical-ml/ClassicalML/run_classical_ml.py \
    --input $inputdir/cts_per_gene_Polyphen2_D.csv \
    --label $labeldir/bbcar_redcap_label_studyid.csv \
    --outfn $outdir/bbcar_genesommut_roc_auc.csv \
    --indexdir $ixdir \
    --scoring roc_auc \
    --savemodel "/projects/b1131/saya/bbcar/models"
#
