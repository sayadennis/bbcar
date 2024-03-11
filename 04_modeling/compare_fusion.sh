#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=4G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="compfusion"
#SBATCH --output=/projects/b1131/saya/bbcar/out/compare_fusion_nested.out

module purge all
module load python-miniconda3/4.12.0
source activate classical-ml

cd /projects/b1131/saya/bbcar/

datadir='/projects/b1131/saya/bbcar/data'
labeldir=$datadir/clinical
outdir=/projects/b1131/saya/bbcar/model_interpretations/breast_cancer_prediction

###################
#### No fusion ####
###################

## Mutation

inputdir=$datadir/02a_mutation/08_feature_matrix/20230423_signature_results
inputfn='sbs_96_original_per_sample.csv'
shortfn=$(echo $fn | cut -d. -f1)
outfn='no_fusion_mutation_orig.csv'
mkdir -p ${outdir}/${shortfn}
python ~/classical-ml/ClassicalML/nested.py \
    ${inputdir}/${inputfn} \
    ${labeldir}/bbcar_label_studyid_from_gatk_filenames.csv \
    ${outdir}/${shortfn}/${outfn} \
    ${SLURM_NTASKS};

## CNV

inputdir=$datadir/02b_cnv/04_cleaned_cnv/all
inputfn='cyto_copy_conf90_all.csv'
shortfn=$(echo $fn | cut -d. -f1)
outfn='no_fusion_cnv_gistic.csv'
mkdir -p ${outdir}/${shortfn}
python ~/classical-ml/ClassicalML/nested.py \
    ${inputdir}/${inputfn} \
    ${labeldir}/bbcar_label_studyid_from_gatk_filenames.csv \
    ${outdir}/${shortfn}/${outfn} \
    ${SLURM_NTASKS};

######################
#### Early fusion ####
######################

inputdir=$datadir/combined_mutation_cnv
inputfn='combined_mut_orig96_cnv_gisticcyto.csv'
shortfn=$(echo $fn | cut -d. -f1)
outfn='early_fusion_orig_with_gisticcyto.csv'
mkdir -p ${outdir}/${shortfn}
python ~/classical-ml/ClassicalML/nested.py \
    ${inputdir}/${inputfn} \
    ${labeldir}/bbcar_label_studyid_from_gatk_filenames.csv \
    ${outdir}/${shortfn}/${outfn} \
    ${SLURM_NTASKS};

#####################
#### Late fusion ####
#####################

input1dir=$datadir/02a_mutation/08_feature_matrix/20230423_signature_results
input1fn='sbs_96_original_per_sample.csv'
input2dir=$datadir/02b_cnv/04_cleaned_cnv/all
input2fn='cyto_copy_conf90_all.csv'
shortfn='late_fusion'
outfn='late_fusion.csv'
mkdir -p ${outdir}/${shortfn}
python ~/classical-ml/ClassicalML/late_fusion.py \
    ${input1dir}/${input1fn} \
    ${input2dir}/${input2fn} \
    ${labeldir}/bbcar_label_studyid_from_gatk_filenames.csv \
    ${outdir}/${shortfn}/${outfn} \
    ${SLURM_NTASKS};

##########################
#### Mid-level fusion ####
##########################

# ## Unsupervised hNMF
# 
# inputdir='/projects/b1131/saya/bbcar/model_interpretations/breast_cancer_prediction/unsupervised_hNMF'
# inputfn='learned_W.csv'
# shortfn=$(echo $fn | cut -d. -f1)
# outfn='mid_fusion_unsupervised.csv'
# mkdir -p ${outdir}/${shortfn}
# python ~/classical-ml/ClassicalML/nested.py \
#     ${inputdir}/${inputfn} \
#     ${labeldir}/bbcar_label_studyid_from_gatk_filenames.csv \
#     ${outdir}/${shortfn}/${outfn} \
#     ${SLURM_NTASKS};
# 

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

