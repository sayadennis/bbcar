#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=3G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="suphNMFml"
#SBATCH --output=/projects/b1131/saya/bbcar/out/supervised_hNMF_W_classicalml.out

module purge all
module load python-miniconda3/4.12.0
source activate classical-ml

cd /projects/b1131/saya/bbcar/

inputdir='/projects/b1131/saya/bbcar/data/combined_mutation_cnv'
labeldir='/projects/b1131/saya/bbcar/data/clinical'
outdir='/projects/b1131/saya/bbcar/model_interpretations'
ixdir='/projects/b1131/saya/bbcar/train_test_splits'

fn='test_learned_W.csv'
shortfn=$(echo $fn | cut -d. -f1)
mkdir -p ${outdir}/${shortfn}
python ~/classical-ml/ClassicalML/run_classical_ml.py \
    --input ${inputdir}/${fn} \
    --label ${labeldir}/bbcar_label_studyid_from_gatk_filenames.csv \
    --outdir ${outdir}/${shortfn}/ \
    --indexdir $ixdir \
    --scoring roc_auc \
    --n_cpu ${SLURM_NTASKS};

