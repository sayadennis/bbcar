#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-gpu
#SBATCH --gres=gpu:a100:1
#SBATCH -t 48:00:00
#SBATCH -n 12
#SBATCH --mem=0
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="scanmap_rerun"
#SBATCH --output=/home/srd6051/bbcar/out/bbcar_run_scanmap_reruns.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

# #### Region-level features #### 

# # Clinical
# python bbcar/scanmap/src/bbcar_run_scanmap.py \
#     --genomic /projects/b1122/saya/scanmap_data/reg_thres.pik \
#     --confound /projects/b1122/saya/scanmap_data/bbcar_clin_intindex.csv \
#     --label /projects/b1122/saya/scanmap_data/bbcar_label_intindex.csv \
#     --outdir /projects/b1042/ClareLab/saya/train_record/regthres_clin \
#     --index /projects/b1122/saya/scanmap_data \
#     --n_iter 4000 \
#     > bbcar/scanmap/results/results_regthres_clin.txt
# #

# # Clinical + driver somatic 
# python bbcar/scanmap/src/bbcar_run_scanmap.py \
#     --genomic /projects/b1122/saya/scanmap_data/reg_thres.pik \
#     --confound /projects/b1122/saya/scanmap_data/bbcar_clin_driversomatic_intindex.csv \
#     --label /projects/b1122/saya/scanmap_data/bbcar_label_intindex.csv \
#     --outdir /projects/b1042/ClareLab/saya/train_record/regthres_clin_driversomatic \
#     --index /projects/b1122/saya/scanmap_data \
#     --n_iter 4000 \
#     > bbcar/scanmap/results/results_regthres_clin_driversomatic.txt
# #

#### Gene-level features #### 

# Clinical
python bbcar/scanmap/src/bbcar_run_scanmap.py \
    --genomic /projects/b1122/saya/scanmap_data/gene_thres_thres010_abs.pik \
    --confound /projects/b1122/saya/scanmap_data/bbcar_clin_intindex.csv \
    --label /projects/b1122/saya/scanmap_data/bbcar_label_intindex.csv \
    --outdir /projects/b1042/ClareLab/saya/train_record/genethres_clin \
    --index /projects/b1122/saya/scanmap_data \
    --n_iter 4000 \
    > bbcar/scanmap/results/results_genethres_clin.txt
#

# Clinical + driver somatic 
python bbcar/scanmap/src/bbcar_run_scanmap.py \
    --genomic /projects/b1122/saya/scanmap_data/gene_thres_thres010_abs.pik \
    --confound /projects/b1122/saya/scanmap_data/bbcar_clin_driversomatic_intindex.csv \
    --label /projects/b1122/saya/scanmap_data/bbcar_label_intindex.csv \
    --outdir /projects/b1042/ClareLab/saya/train_record/genethres_clin_driversomatic \
    --index /projects/b1122/saya/scanmap_data \
    --n_iter 4000 \
    > bbcar/scanmap/results/results_genethres_clin_driversomatic.txt
#

