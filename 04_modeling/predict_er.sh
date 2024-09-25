#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=4G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="predER"
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/predict_er.out

module purge all
module load python-miniconda3/4.12.0
source activate classical-ml

cd /projects/b1131/saya/new_bbcar/

datadir='/projects/b1131/saya/new_bbcar/data'
labeldir='/projects/b1131/saya/new_bbcar'
outdir='/projects/b1131/saya/new_bbcar/model_interpretations/breast_cancer_prediction/predict_er'

mkdir -p $outdir

###################################
#### Including controls as HR- ####
###################################

## Mutation

inputdir=$datadir/02a_mutation/08_feature_matrix
inputfn='raw_SBS96_features.csv'
outfn='mutation_with_cont.csv'
python ~/classical-ml/ClassicalML/nested.py \
    ${inputdir}/${inputfn} \
    ${labeldir}/label_er_with_cont.csv \
    ${outdir}/${outfn} \
    ${SLURM_NTASKS};

## CNV

inputdir=$datadir/02b_cnv/10_cleaned_cnv
inputfn='cyto_thres_aber.csv'
outfn='cnv_with_cont.csv'
python ~/classical-ml/ClassicalML/nested.py \
    ${inputdir}/${inputfn} \
    ${labeldir}/label_er_with_cont.csv \
    ${outdir}/${outfn} \
    ${SLURM_NTASKS};

############################
#### Excluding controls ####
############################

## Mutation

inputdir=$datadir/02a_mutation/08_feature_matrix
inputfn='raw_SBS96_features.csv'
outfn='mutation_without_cont.csv'
python ~/classical-ml/ClassicalML/nested.py \
    ${inputdir}/${inputfn} \
    ${labeldir}/label_er_without_cont.csv \
    ${outdir}/${outfn} \
    ${SLURM_NTASKS};

## CNV

inputdir=$datadir/02b_cnv/10_cleaned_cnv
inputfn='cyto_thres_aber.csv'
outfn='cnv_without_cont.csv'
python ~/classical-ml/ClassicalML/nested.py \
    ${inputdir}/${inputfn} \
    ${labeldir}/label_er_without_cont.csv \
    ${outdir}/${outfn} \
    ${SLURM_NTASKS};

