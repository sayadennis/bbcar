#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-gpu
#SBATCH --gres=gpu:a100:1
#SBATCH -t 48:00:00
#SBATCH -n 12
#SBATCH --mem=0
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="scanmap_siamese_nocf"
#SBATCH --output=bbcar/out/bbcar_run_scanmap_siamese_nocf.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

datadir="/projects/b1122/saya/06_modified_data"
labeldir="/projects/b1122/saya"
train_record_dir="/projects/b1042/ClareLab/saya/train_record_siamese_nocf"
indexdir="/projects/b1122/saya/indices"
resultsdir="bbcar/model_performance/results_scanmap_siamese_nocf"

#### Region-level features #### 

python bbcar/src/modeling/exploratory_models/bbcar_run_scanmap_siamese_nocf.py \
    --genomic $datadir/reg_thres_conf90_intindex.csv \
    --label $labeldir/bbcar_label_intindex.csv \
    --outdir $train_record_dir \
    --index $indexdir \
    --n_iter 4000 \
    > $resultsdir/results_scanmap_siamese_regthres_nocf.txt
#

#### Gene-level features #### 

python bbcar/src/modeling/exploratory_models/bbcar_run_scanmap_siamese_nocf.py \
    --genomic $datadir/gene_thres_conf90_intindex.csv \
    --label $labeldir/bbcar_label_intindex.csv \
    --outdir $train_record_dir \
    --index $indexdir \
    --n_iter 4000 \
    > $resultsdir/results_scanmap_siamese_genethres_nocf.txt
#

