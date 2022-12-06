#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=1G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="cnsigml"
#SBATCH --output=/projects/b1131/saya/bbcar/out/cn_signature_classicalml.out

module purge all
module load python-miniconda3/4.12.0
source activate classical-ml

cd /projects/b1131/saya/bbcar/

inputdir='/projects/b1131/saya/bbcar/data/02b_cnv/signatures'
labeldir='/projects/b1131/saya/bbcar/data/clinical'
outdir='/projects/b1131/saya/bbcar/model_interpretations'
ixdir='/projects/b1131/saya/bbcar/train_test_splits'

## Non-NMF ML training
for fn in $(ls /projects/b1131/saya/bbcar/data/02b_cnv/signatures/); do
    shortfn=$(echo $fn | cut -d. -f1)
    mkdir -p ${outdir}/cnv_sig_nonNMF/${shortfn}
    python ~/classical-ml/ClassicalML/run_classical_ml.py \
        --input ${inputdir}/${fn} \
        --label ${labeldir}/bbcar_label_studyid_from_gatk_filenames.csv \
        --outdir ${outdir}/cnv_sig_nonNMF/${shortfn}/ \
        --indexdir $ixdir \
        --scoring roc_auc;
done

## Non-NMF ML training
for fn in $(ls /projects/b1131/saya/bbcar/data/02b_cnv/signatures/); do
    shortfn=$(echo $fn | cut -d. -f1)
    mkdir -p ${outdir}/cnv_sig_NMF/${shortfn}
    python ~/classical-ml/ClassicalML/run_classical_ml.py \
        --input ${inputdir}/${fn} \
        --label ${labeldir}/bbcar_label_studyid_from_gatk_filenames.csv \
        --outdir ${outdir}/cnv_sig_NMF/${shortfn}/ \
        --indexdir $ixdir \
        --scoring roc_auc \
        --nmf 500;
done
