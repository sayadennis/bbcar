#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -n 12
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="bbcarcnv_ccsep"
#SBATCH --output=bbcar/out/bbcar_cnv_cytoband_gene_ccsep.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

inputdir='/projects/b1122/saya/06_modified_data'
labeldir='/projects/b1122/saya'
outdir='/home/srd6051/bbcar/model_performance/results_classicalml'
ixdir='/projects/b1122/saya/indices'

## Modeling cytoband-level data

# tune using accuracy 
python classical_ml/ClassicalML/run_classical_ml.py \
    --input $inputdir/cyto_copy_conf90_ccsep_studyindex.csv \
    --label $labeldir/bbcar_label_studyindex.csv \
    --outfn $outdir/bbcar_cytocopy_ccsep_accuracy.csv \
    --indexdir $ixdir \
    --scoring accuracy

# tune using ROC-AUC
python classical_ml/ClassicalML/run_classical_ml.py \
    --input $inputdir/cyto_copy_conf90_ccsep_studyindex.csv \
    --label $labeldir/bbcar_label_studyindex.csv \
    --outfn $outdir/bbcar_cytocopy_ccsep_roc_auc.csv \
    --indexdir $ixdir \
    --scoring roc_auc

# tune using precision 
python classical_ml/ClassicalML/run_classical_ml.py \
    --input $inputdir/cyto_copy_conf90_ccsep_studyindex.csv \
    --label $labeldir/bbcar_label_studyindex.csv \
    --outfn $outdir/bbcar_cytocopy_ccsep_precision.csv \
    --indexdir $ixdir \
    --scoring precision

# tune using f1 
python classical_ml/ClassicalML/run_classical_ml.py \
    --input $inputdir/cyto_copy_conf90_ccsep_studyindex.csv \
    --label $labeldir/bbcar_label_studyindex.csv \
    --outfn $outdir/bbcar_cytocopy_ccsep_f1.csv \
    --indexdir $ixdir \
    --scoring f1
#

## Modeling gene-level data 

# tune using precision 
python classical_ml/ClassicalML/run_classical_ml.py \
    --input $inputdir/gene_copy_conf90_ccsep_studyindex.csv \
    --label $labeldir/bbcar_label_studyindex.csv \
    --outfn $outdir/bbcar_genecopy_ccsep_precision.csv \
    --indexdir $ixdir \
    --scoring precision
#
