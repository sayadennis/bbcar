#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-gpu
#SBATCH --gres=gpu:a100:1
#SBATCH -t 48:00:00
#SBATCH -n 12
#SBATCH --mem=0
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="regclin_siamese"
#SBATCH --output=bbcar/out/bbcar_run_scanmap_siamese_regthres_clin.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv

datadir="/projects/b1122/saya/06_modified_data"
cfdir="/projects/b1122/saya/bbcar_non_cnv_features"
labeldir="/projects/b1122/saya"
train_record_dir="/projects/b1042/ClareLab/saya/train_record_siamese"
indexdir="/projects/b1122/saya/indices"
resultsdir="bbcar/model_performance/results_scanmap_siamese"

#### Region-level features #### 

# Clinical
python bbcar/src/modeling/scanmap/src/bbcar_run_scanmap_siamese.py \
    --genomic $datadir/reg_thres_conf90_intindex.csv \
    --confound $cfdir/bbcar_clin_intindex.csv \
    --label $labeldir/bbcar_label_intindex.csv \
    --outdir $train_record_dir/regthres_clin \
    --index $indexdir \
    --n_iter 4000 \
    > $resultsdir/results_siamese_regthres_clin.txt
#

# Clinical + mutational signature 
python bbcar/src/modeling/scanmap/src/bbcar_run_scanmap_siamese.py \
    --genomic $datadir/reg_thres_conf90_intindex.csv \
    --confound $cfdir/bbcar_clin_mut_intindex.csv \
    --label $labeldir/bbcar_label_intindex.csv \
    --outdir $train_record_dir/regthres_clin_mut \
    --index $indexdir \
    > $resultsdir/results_regthres_clin_mut.txt
#

# Clinical + PRS 
python bbcar/src/modeling/scanmap/src/bbcar_run_scanmap_siamese.py \
    --genomic $datadir/reg_thres_conf90_intindex.csv \
    --confound $cfdir/bbcar_clin_prs_intindex.csv \
    --label $labeldir/bbcar_label_intindex.csv \
    --outdir $train_record_dir/regthres_clin_prs \
    --index $indexdir \
    > $resultsdir/results_regthres_clin_prs.txt
#

# Clinical + driver somatic 
python bbcar/src/modeling/scanmap/src/bbcar_run_scanmap_siamese.py \
    --genomic $datadir/reg_thres_conf90_intindex.csv \
    --confound $cfdir/bbcar_clin_driversomatic_intindex.csv \
    --label $labeldir/bbcar_label_intindex.csv \
    --outdir $train_record_dir/regthres_clin_driversomatic \
    --index $indexdir \
    --n_iter 4000 \
    > $resultsdir/results_regthres_clin_driversomatic.txt
#

# Clinical + mutational signature + PRS
python bbcar/src/modeling/scanmap/src/bbcar_run_scanmap_siamese.py \
    --genomic $datadir/reg_thres_conf90_intindex.csv \
    --confound $cfdir/bbcar_clin_mut_prs_intindex.csv \
    --label $labeldir/bbcar_label_intindex.csv \
    --outdir $train_record_dir/regthres_clin_mut_prs \
    --index $indexdir \
    > $resultsdir/results_regthres_clin_mut_prs.txt
#

# Clinical + mutational signature + driver somatic mutation
python bbcar/src/modeling/scanmap/src/bbcar_run_scanmap_siamese.py \
    --genomic $datadir/reg_thres_conf90_intindex.csv \
    --confound $cfdir/bbcar_clin_mut_driversomatic_intindex.csv \
    --label $labeldir/bbcar_label_intindex.csv \
    --outdir $train_record_dir/regthres_clin_mut_driversomatic \
    --index $indexdir \
    > $resultsdir/results_regthres_clin_mut_driversomatic.txt
#

# Clinical + PRS + driver somatic mutation 
python bbcar/src/modeling/scanmap/src/bbcar_run_scanmap_siamese.py \
    --genomic $datadir/reg_thres_conf90_intindex.csv \
    --confound $cfdir/bbcar_clin_prs_driversomatic_intindex.csv \
    --label $labeldir/bbcar_label_intindex.csv \
    --outdir $train_record_dir/regthres_clin_prs_driversomatic \
    --index $indexdir \
    > $resultsdir/results_regthres_clin_prs_driversomatic.txt
#
