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
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/compare_gistic_outs.out

module purge all
module load python-miniconda3/4.12.0
source activate classical-ml

cd /projects/b1131/saya/new_bbcar/

inputdir='/projects/b1131/saya/new_bbcar/data/02b_cnv/11_batcheffect_removed'
labeldir='/projects/b1131/saya/new_bbcar'
outdir='/projects/b1131/saya/new_bbcar/model_interpretations/breast_cancer_prediction/compare_gistic'

mkdir -p $outdir

for inputfn in "cyto_copy.csv" "cyto_thres_aber.csv" "reg_copy.csv" "reg_thres.csv"; do
    python ~/classical-ml/ClassicalML/nested.py \
        ${inputdir}/${inputfn} \
        ${labeldir}/label_all.csv \
        ${outdir}/results_${inputfn} \
        ${SLURM_NTASKS};
done

