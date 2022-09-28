#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-gpu
#SBATCH --gres=gpu:a100:1
#SBATCH -t 48:00:00
#SBATCH -n 12
#SBATCH --mem=0
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="regmut_siamese"
#SBATCH --output=bbcar/out/bbcar_run_scanmap_siamese_regthres_mut.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

datadir="/projects/b1122/saya/06_modified_data"
cfdir="/projects/b1122/saya/bbcar_non_cnv_features"
labeldir="/projects/b1122/saya"
train_record_dir="/projects/b1042/ClareLab/saya/train_record_siamese"
indexdir="/projects/b1122/saya/indices"
resultsdir="bbcar/model_performance/results_scanmap_siamese"

#### Region-level features #### 

# Only use mutational signature as confounding variables 
python bbcar/src/modeling/scanmap/src/bbcar_run_scanmap_siamese.py \
    --genomic $datadir/reg_thres_conf90_intindex.csv \
    --confound $cfdir/bbcar_mut_intindex.csv \
    --label $labeldir/bbcar_label_intindex.csv \
    --outdir $train_record_dir/regthres_mut \
    --index $indexdir \
    --n_iter 4000 \
    > $resultsdir/results_siamese_regthres_mut.txt
#
