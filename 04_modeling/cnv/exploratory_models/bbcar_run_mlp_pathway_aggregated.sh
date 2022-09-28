#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-gpu
#SBATCH --gres=gpu:a100:1
#SBATCH -t 48:00:00
#SBATCH -n 12
#SBATCH --mem=0
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="bbcarmlp_pathway"
#SBATCH --output=bbcar/out/bbcar_run_mlp_pathway_aggregated.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

datadir="/projects/b1122/saya/08_pathway_aggregated_cnv"
labeldir="/projects/b1122/saya"
train_record_dir="/projects/b1042/ClareLab/saya/train_record_mlp/pathway_aggregated"
indexdir="/projects/b1122/saya/indices"
resultsdir="bbcar/model_performance/results_mlp"

# Averaged 
python bbcar/src/modeling/exploratory_models/bbcar_run_mlp.py \
    --genomic $datadir/pathway_avg_cnv_intindex.csv \
    --label $labeldir/bbcar_label_intindex.csv \
    --outdir $train_record_dir/avg \
    --index $indexdir \
    --n_iter 2000 \
    > $resultsdir/results_pathway_avg.txt
#

# Max
python bbcar/src/modeling/exploratory_models/bbcar_run_mlp.py \
    --genomic $datadir/pathway_max_cnv_intindex.csv \
    --label $labeldir/bbcar_label_intindex.csv \
    --outdir $train_record_dir/max \
    --index $indexdir \
    --n_iter 2000 \
    > $resultsdir/results_pathway_max.txt
#

# Median
python bbcar/src/modeling/exploratory_models/bbcar_run_mlp.py \
    --genomic $datadir/pathway_med_cnv_intindex.csv \
    --label $labeldir/bbcar_label_intindex.csv \
    --outdir $train_record_dir/med \
    --index $indexdir \
    --n_iter 2000 \
    > $resultsdir/results_pathway_med.txt
#

# Sum
python bbcar/src/modeling/exploratory_models/bbcar_run_mlp.py \
    --genomic $datadir/pathway_sum_cnv_intindex.csv \
    --label $labeldir/bbcar_label_intindex.csv \
    --outdir $train_record_dir/sum \
    --index $indexdir \
    --n_iter 2000 \
    > $resultsdir/results_pathway_sum.txt
#
