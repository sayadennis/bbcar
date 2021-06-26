#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-gpu
#SBATCH --gres=gpu:a100:1
#SBATCH -t 48:00:00
#SBATCH -n 12
#SBATCH --mem=0
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="bbcarcontmlp"
#SBATCH --output=bbcar/out/bbcar_run_mlp_contrastive.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

datadir="/projects/b1122/saya/06_modified_data"
labeldir="/projects/b1122/saya"
train_record_dir="/projects/b1042/ClareLab/saya/train_record_mlp_contrastive"
indexdir="/projects/b1122/saya/indices"
resultsdir="bbcar/model_performance/results_mlp_contrastive"

#### Region-level features #### 

# Region-level threshold data
python bbcar/src/modeling/exploratory_models/bbcar_run_mlp.py \
    --genomic $datadir/reg_thres_conf90_intindex.csv \
    --label $labeldir/bbcar_label_intindex.csv \
    --outdir $train_record_dir \
    --index $indexdir \
    --contrastive euclidean \
    --n_iter 2000 \
    > $resultsdir/results_regthres.txt
#

# Region-level copy number data
python bbcar/src/modeling/exploratory_models/bbcar_run_mlp.py \
    --genomic $datadir/reg_copy_conf90_intindex.csv \
    --label $labeldir/bbcar_label_intindex.csv \
    --outdir $train_record_dir \
    --index $indexdir \
    --contrastive euclidean \
    --n_iter 2000 \
    > $resultsdir/results_regcopy.txt
#

#### Gene-level features #### 

# Gene-level threshold data
python bbcar/src/modeling/exploratory_models/bbcar_run_mlp.py \
    --genomic $datadir/gene_thres_conf90_intindex.csv \
    --label $labeldir/bbcar_label_intindex.csv \
    --outdir $train_record_dir \
    --index $indexdir \
    --contrastive euclidean \
    --n_iter 2000 \
    > $resultsdir/results_genethres.txt
#

# Gene-level copy number data
python bbcar/src/modeling/exploratory_models/bbcar_run_mlp.py \
    --genomic $datadir/gene_copy_conf90_intindex.csv \
    --label $labeldir/bbcar_label_intindex.csv \
    --outdir $train_record_dir \
    --index $indexdir \
    --contrastive euclidean \
    --n_iter 2000 \
    > $resultsdir/results_genecopy.txt
#
