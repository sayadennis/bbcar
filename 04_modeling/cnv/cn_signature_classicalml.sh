#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=3G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="cnsigml"
#SBATCH --output=/projects/b1131/saya/bbcar/out/cn_inhouse_signature_classicalml.out

module purge all
module load python-miniconda3/4.12.0
source activate classical-ml

cd /projects/b1131/saya/bbcar/

inputdir='/projects/b1131/saya/bbcar/data/02b_cnv/inhouse_signatures'
labeldir='/projects/b1131/saya/bbcar/data/clinical'
outdir='/projects/b1131/saya/bbcar/model_interpretations/breast_cancer_prediction'
ixdir='/projects/b1131/saya/bbcar/train_test_splits'

## ML modeling on original CN features
fn='inhouse_sig_batcheffect_rm_combat.csv'
shortfn=$(echo $fn | cut -d. -f1)
mkdir -p ${outdir}/cnv_sig_nonNMF/${shortfn}
python ~/classical-ml/ClassicalML/run_classical_ml.py \
    --input ${inputdir}/${fn} \
    --label ${labeldir}/bbcar_label_studyid_from_gatk_filenames.csv \
    --outdir ${outdir}/cnv_sig_nonNMF/${shortfn}/ \
    --indexdir $ixdir \
    --scoring roc_auc \
    --n_cpu ${SLURM_NTASKS};
#

## ML modeling on CN features that have been converted to signatures via NMF
fn='inhouse_cnv_sig_per_sample.csv'
shortfn=$(echo $fn | cut -d. -f1)
mkdir -p ${outdir}/cnv_sig_nonNMF/${shortfn}
python ~/classical-ml/ClassicalML/run_classical_ml.py \
    --input ${inputdir}/${fn} \
    --label ${labeldir}/bbcar_label_studyid_from_gatk_filenames.csv \
    --outdir ${outdir}/cnv_sig_nonNMF/${shortfn}/ \
    --indexdir $ixdir \
    --scoring roc_auc \
    --n_cpu ${SLURM_NTASKS};
#

### NMF ML training
#fn='inhouse_sig_batcheffect_rm_combat.csv'
#shortfn=$(echo $fn | cut -d. -f1)
#mkdir -p ${outdir}/cnv_sig_NMF/${shortfn}
#python ~/classical-ml/ClassicalML/run_classical_ml.py \
#    --input ${inputdir}/${fn} \
#    --label ${labeldir}/bbcar_label_studyid_from_gatk_filenames.csv \
#    --outdir ${outdir}/cnv_sig_NMF/${shortfn}/ \
#    --indexdir $ixdir \
#    --scoring roc_auc \
#    --n_cpu ${SLURM_NTASKS} \
#    --nmf 500;

