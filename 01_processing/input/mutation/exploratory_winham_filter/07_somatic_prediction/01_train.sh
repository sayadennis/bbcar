#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicslong
#SBATCH -t 240:00:00
#SBATCH -N 1
#SBATCH -n 48
#SBATCH --mem=96G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="bbcarsom_ml"
#SBATCH --output=/projects/b1131/saya/bbcar/out/bbcar_somatic_classicalml.out

module purge all
module load python-miniconda3/4.12.0
source activate classical-ml

cd /projects/b1131/saya/bbcar

pon_source='bbcar'

for filter_type in liberal classical strict; do
    inputdir=/projects/b1131/saya/bbcar/exploratory_winham_filter/${filter_type}/04_ml_features
    labeldir=/projects/b1131/saya/bbcar/exploratory_winham_filter/${filter_type}/04_ml_features
    ixdir=/projects/b1131/saya/bbcar/exploratory_winham_filter/${filter_type}/04_ml_features/somatic_pred_ix/${pon_source}
    outdir=/projects/b1131/saya/bbcar/model_interpretations/exploratory_winham_filter/${filter_type}
    # modeldir=/projects/b1131/saya/bbcar/models/exploratory_winham_filter/${filter_type}

    mkdir -p $outdir
    # mkdir -p $modeldir

    python ${HOME}/classical-ml/ClassicalML/run_classical_ml.py \
        --input $inputdir/input_matched_${pon_source}.csv \
        --label $labeldir/target_matched_${pon_source}.csv \
        --outdir $outdir \
        --indexdir $ixdir \
        --scoring roc_auc \
        --n_cpu ${SLURM_NTASKS}
done
