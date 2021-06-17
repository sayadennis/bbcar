#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-gpu
#SBATCH --gres=gpu:a100:1
#SBATCH -t 48:00:00
#SBATCH -n 12
#SBATCH --mem=0
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="bbcarmlp"
#SBATCH --output=bbcar/out/bbcar_run_mlp.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

#### Region-level features #### 

# Region-level threshold data
python bbcar/src/modeling/exploratory_models/bbcar_run_mlp.py \
    --genomic /projects/b1122/saya/06_modified_data/reg_thres_conf90_intindex.csv \
    --label /projects/b1122/saya/bbcar_label_intindex.csv \
    --outdir /projects/b1042/ClareLab/saya/train_record_mlp \
    --index /projects/b1122/saya/indices \
    --n_iter 2000 \
    > bbcar/model_performance/results_mlp/results_regthres.txt
#

# Region-level copy number data
python bbcar/src/modeling/exploratory_models/bbcar_run_mlp.py \
    --genomic /projects/b1122/saya/06_modified_data/reg_copy_conf90_intindex.csv \
    --label /projects/b1122/saya/bbcar_label_intindex.csv \
    --outdir /projects/b1042/ClareLab/saya/train_record_mlp \
    --index /projects/b1122/saya/indices \
    --n_iter 2000 \
    > bbcar/model_performance/results_mlp/results_regcopy.txt
#

#### Gene-level features #### 

# Gene-level threshold data
python bbcar/src/modeling/exploratory_models/bbcar_run_mlp.py \
    --genomic /projects/b1122/saya/06_modified_data/gene_thres_conf90_intindex.csv \
    --label /projects/b1122/saya/bbcar_label_intindex.csv \
    --outdir /projects/b1042/ClareLab/saya/train_record_mlp \
    --index /projects/b1122/saya/indices \
    --n_iter 2000 \
    > bbcar/model_performance/results_mlp/results_genethres.txt
#

# Gene-level copy number data
python bbcar/src/modeling/exploratory_models/bbcar_run_mlp.py \
    --genomic /projects/b1122/saya/06_modified_data/gene_copy_conf90_intindex.csv \
    --label /projects/b1122/saya/bbcar_label_intindex.csv \
    --outdir /projects/b1042/ClareLab/saya/train_record_mlp \
    --index /projects/b1122/saya/indices \
    --n_iter 2000 \
    > bbcar/model_performance/results_mlp/results_genecopy.txt
#
