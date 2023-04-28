#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=16G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="sigpro_matched"
#SBATCH --output=/projects/b1131/saya/bbcar/out/run_SigProfilerExtractor_matched.out

module purge all
module load python-miniconda3/4.12.0
source activate SigProfilerExtractor

cd /projects/b1131/saya/bbcar/

pon_source=bbcar
for filter_type in liberal classical strict; do
    mkdir -p /projects/b1131/saya/bbcar/exploratory_winham_filter/${filter_type}/08_feature_matrix/signature_results
    python ~/bbcar/repo/01_processing/input/mutation/exploratory_winham_filter/08_feature_matrix/run_SigProfilerExtractor.py $filter_type
done
