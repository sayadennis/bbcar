#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=4G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="testnmfcv"
#SBATCH --output=/projects/b1131/saya/bbcar/out/test_nmf_nested_cv.out

module purge all
module load python-miniconda3/4.12.0
source activate classical-ml

cd /projects/b1131/saya/bbcar/

datadir='/projects/b1131/saya/bbcar/data'
labeldir=$datadir/clinical
outdir=/projects/b1131/saya/bbcar/model_interpretations/breast_cancer_prediction

##########################
#### Mid-level fusion ####
##########################

## Supervised hNMF
input1dir=$datadir/02a_mutation/08_feature_matrix/20230423_signature_results
input1fn='sbs_96_original_per_sample.csv'
input2dir=$datadir/02b_cnv/inhouse_signatures
input2fn='inhouse_cn_features_batcheffect_rm_combat.csv'
shortfn='mid_fusion_supervised'
outfn='mid_fusion_supervised.csv'
mkdir -p ${outdir}/${shortfn}
python ~/classical-ml/ClassicalML/nmf_nested.py \
    ${input1dir}/${input1fn} \
    ${input2dir}/${input2fn} \
    ${labeldir}/bbcar_label_studyid_from_gatk_filenames.csv \
    ${outdir}/${shortfn}/${outfn} \
    ${SLURM_NTASKS};

