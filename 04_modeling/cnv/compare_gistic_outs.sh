#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=4G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=compgistic
#SBATCH --output=/projects/b1131/saya/bbcar/out/compare_gistic_outs.out

module purge all
module load python-miniconda3/4.12.0
source activate classical-ml

cd /projects/b1131/saya/bbcar/

datadir='/projects/b1131/saya/bbcar/data'
labeldir=$datadir/clinical
inputdir=$datadir/02b_cnv/10_cleaned_cnv
outdir=/projects/b1131/saya/bbcar/model_interpretations/breast_cancer_prediction/compare_gistic

mkdir -p $outdir

for inputfn in "cyto_copy_conf90.csv" "cyto_thres_aber_conf90.csv" "reg_copy_conf90.csv" "reg_thres_conf90.csv"; do
    python ~/classical-ml/ClassicalML/nested.py \
        ${inputdir}/${inputfn} \
        ${labeldir}/bbcar_label_studyid_from_gatk_filenames.csv \
        ${outdir}/results_${inputfn} \
        ${SLURM_NTASKS};
done
