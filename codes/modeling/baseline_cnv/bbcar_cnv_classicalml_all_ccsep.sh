#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="bbcar_cnv_classicalml_all_ccsep"
#SBATCH --output=bbcar/out/bbcar_cnv_classicalml_all_ccsep.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

inputdir='/projects/b1122/saya/04_cleaned_cnv/all_ccsep'
labeldir='/projects/b1122/saya'
outdir='/home/srd6051/bbcar/model_performance/results_classicalml/all_ccsep'
ixdir='/projects/b1122/saya/indices'

#############################
#### Gene-level features ####
#############################

## Gene-level copy number
python classical_ml/ClassicalML/run_classical_ml.py \
    --input $inputdir/gene_copy_conf90_ccsep.csv \
    --label $labeldir/bbcar_label_studyid.csv \
    --outfn $outdir/bbcar_genecopy_roc_auc_ovo.csv \
    --indexdir $ixdir \
    --scoring roc_auc_ovo
#

## Gene-level threshold values
python classical_ml/ClassicalML/run_classical_ml.py \
    --input $inputdir/gene_thres_conf90_ccsep.csv \
    --label $labeldir/bbcar_label_studyid.csv \
    --outfn $outdir/bbcar_genethres_roc_auc_ovo.csv \
    --indexdir $ixdir \
    --scoring roc_auc_ovo
#

###############################
#### Region-level features ####
###############################

## Region-level copy number
python classical_ml/ClassicalML/run_classical_ml.py \
    --input $inputdir/reg_copy_unique.csv \
    --label $labeldir/bbcar_label_studyid.csv \
    --outfn $outdir/bbcar_regcopy_roc_auc_ovo.csv \
    --indexdir $ixdir \
    --scoring roc_auc_ovo

# --mem=0

#################################
#### Cytoband-level features ####
#################################

python classical_ml/ClassicalML/run_classical_ml.py \
    --input $inputdir/cyto_copy_conf90_ccsep.csv \
    --label $labeldir/bbcar_label_studyid.csv \
    --outfn $outdir/bbcar_cytocopy_roc_auc_ovo.csv \
    --indexdir $ixdir \
    --scoring roc_auc_ovo
#
