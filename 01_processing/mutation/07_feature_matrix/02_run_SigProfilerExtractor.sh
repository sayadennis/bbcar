#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=24G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="sigpro"
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/run_SigProfilerExtractor_rm10k.out

module purge all
module load python-miniconda3/4.12.0
source activate SigProfilerExtractor

mutation_dir='/projects/b1131/saya/new_bbcar/data/02a_mutation'

cd "${mutation_dir}/08_feature_matrix/"

python ~/bbcar/repo/01_processing/mutation/07_feature_matrix/run_SigProfilerExtractor.py \
    "${mutation_dir}/07_predicted_somatic/vcfs_rm10k" \
    "${mutation_dir}/08_feature_matrix/signature_results_rm10k" \
    ${SLURM_NTASKS}

