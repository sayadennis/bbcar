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
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/compare_fusion.out

module purge all
module load python-miniconda3/4.12.0
source activate classical-ml

cd /projects/b1131/saya/new_bbcar/

datadir='/projects/b1131/saya/new_bbcar/data'
labeldir='/projects/b1131/saya/new_bbcar'
outdir='/projects/b1131/saya/new_bbcar/model_interpretations/breast_cancer_prediction'

# ###################
# #### No fusion ####
# ###################
# 
# ## Mutation
# 
# inputdir=$datadir/02a_mutation/08_feature_matrix
# inputfn='raw_SBS96_features.csv'
# outfn='no_fusion_mutation_orig.csv'
# python ~/classical-ml/ClassicalML/nested.py \
#     ${inputdir}/${inputfn} \
#     ${labeldir}/label_all.csv \
#     ${outdir}/${outfn} \
#     ${SLURM_NTASKS};
# 
# ## CNV
# 
# inputdir=$datadir/02b_cnv/10_cleaned_cnv
# inputfn='cyto_thres_aber.csv'
# outfn='no_fusion_cnv_cytothres.csv'
# python ~/classical-ml/ClassicalML/nested.py \
#     ${inputdir}/${inputfn} \
#     ${labeldir}/label_all.csv \
#     ${outdir}/${outfn} \
#     ${SLURM_NTASKS};
# 
# ######################
# #### Early fusion ####
# ######################
# 
# inputdir=$datadir
# inputfn='combined_mutation_cnv_cytothres.csv'
# outfn='early_fusion_cnv_cytothres.csv'
# python ~/classical-ml/ClassicalML/nested.py \
#     ${inputdir}/${inputfn} \
#     ${labeldir}/label_all.csv \
#     ${outdir}/${outfn} \
#     ${SLURM_NTASKS};
# 
# #####################
# #### Late fusion ####
# #####################
# 
# input1dir=$datadir/02a_mutation/08_feature_matrix
# input1fn='raw_SBS96_features.csv'
# input2dir=$datadir/02b_cnv/10_cleaned_cnv
# input2fn='cyto_thres_aber.csv'
# outfn='late_fusion_cnv_cytothres.csv'
# python ~/classical-ml/ClassicalML/late_fusion.py \
#     ${input1dir}/${input1fn} \
#     ${input2dir}/${input2fn} \
#     ${labeldir}/label_all.csv \
#     ${outdir}/${outfn} \
#     ${SLURM_NTASKS};

##########################
#### Mid-level fusion ####
##########################

# ## Unsupervised hNMF
# 
# inputdir='/projects/b1131/saya/bbcar/model_interpretations/breast_cancer_prediction/unsupervised_hNMF'
# inputfn='learned_W.csv'
# outfn='mid_fusion_unsupervised.csv'
# python ~/classical-ml/ClassicalML/nested.py \
#     ${inputdir}/${inputfn} \
#     ${labeldir}/label_all.csv \
#     ${outdir}/${outfn} \
#     ${SLURM_NTASKS};
# 

## Supervised hNMF
input1dir=$datadir/02a_mutation/08_feature_matrix
input1fn='raw_SBS96_features.csv'
input2dir=$datadir/02b_cnv/10_cleaned_cnv
input2fn='cyto_thres_aber.csv'
outfn='mid_fusion_supervised_cnv_cytothres.csv'
python ~/classical-ml/ClassicalML/nmf_nested.py \
    ${input1dir}/${input1fn} \
    ${input2dir}/${input2fn} \
    ${labeldir}/label_all.csv \
    ${outdir}/${outfn} \
    ${SLURM_NTASKS};

