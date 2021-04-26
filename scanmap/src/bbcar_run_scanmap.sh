#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-gpu
#SBATCH --gres=gpu:a100:1
#SBATCH -t 48:00:00
#SBATCH -n 1
#SBATCH --mail-user=<email>
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="scanmap_bbcar"
#SBATCH --output=bbcar/out/210423_scanmap_bbcar.out

. ~/anaconda3/etc/profile.d/conda.sh
conda activate bbcarenv
# --mem=225G ?

#### Region-level features #### 

# Clinical
python bbcar/scanmap/src/bbcar_run_scanmap.py \
    --genomic /projects/b1122/saya/scanmap_data/reg_thres.pik \
    --confound /projects/b1122/saya/scanmap_data/pt_demo_intindex.csv \
    --label /projects/b1122/saya/scanmap_data/bbcar_label_intindex.csv \
    --outdir /projects/b1122/saya/scanmap_data/out \
    --index /projects/b1122/saya/scanmap_data \
    > bbcar/scanmap/results/results_regthres_clin.txt
#

# Clinical + mutational signature 
python bbcar/scanmap/src/bbcar_run_scanmap.py \
    --genomic /projects/b1122/saya/scanmap_data/reg_thres.pik \
    --confound /projects/b1122/saya/scanmap_data/pt_demo_clinmut_intindex.csv \
    --label /projects/b1122/saya/scanmap_data/bbcar_label_intindex.csv \
    --outdir /projects/b1122/saya/scanmap_data/out \
    --index /projects/b1122/saya/scanmap_data \
    > bbcar/scanmap/results/results_regthres_clin_mut.txt
#

# Clinical + mutational signature + somatic mutation in driver genes 
python bbcar/scanmap/src/bbcar_run_scanmap.py \
    --genomic /projects/b1122/saya/scanmap_data/reg_thres.pik \
    --confound /projects/b1122/saya/scanmap_data/pt_demo_clin_mut_driversomatic_intindex.csv \
    --label /projects/b1122/saya/scanmap_data/bbcar_label_intindex.csv \
    --outdir /projects/b1122/saya/scanmap_data/out \
    --index /projects/b1122/saya/scanmap_data \
    > bbcar/scanmap/results/results_regthres_clin_mut_driversomatic.txt
#

# Clinical + PRS SNPs (removed ones appearing in <=3 patients) 
python bbcar/scanmap/src/bbcar_run_scanmap.py \
    --genomic /projects/b1122/saya/scanmap_data/reg_thres.pik \
    --confound /projects/b1122/saya/scanmap_data/clin_prssnp_rm3.csv \
    --label /projects/b1122/saya/scanmap_data/bbcar_label_intindex.csv \
    --outdir /projects/b1122/saya/scanmap_data/out \
    --index /projects/b1122/saya/scanmap_data \
    > bbcar/scanmap/results/results_regthres_clin_prs.txt
#

# Clinical + mutational signature + PRS SNPs 
python bbcar/scanmap/src/bbcar_run_scanmap.py \
    --genomic /projects/b1122/saya/scanmap_data/reg_thres.pik \
    --confound /projects/b1122/saya/scanmap_data/clinmut_prssnp_rm3.csv \
    --label /projects/b1122/saya/scanmap_data/bbcar_label_intindex.csv \
    --outdir /projects/b1122/saya/scanmap_data/out \
    --index /projects/b1122/saya/scanmap_data \
    > bbcar/scanmap/results/results_regthres_clin_mut_prs.txt
#

#### Gene-level features #### 

# Clinical
python bbcar/scanmap/src/bbcar_run_scanmap.py \
    --genomic /projects/b1122/saya/scanmap_data/gene_thres_thres010_abs.pik \
    --confound /projects/b1122/saya/scanmap_data/pt_demo_intindex.csv \
    --label /projects/b1122/saya/scanmap_data/bbcar_label_intindex.csv \
    --outdir /projects/b1122/saya/scanmap_data/out \
    --index /projects/b1122/saya/scanmap_data \
    > bbcar/scanmap/results/results_genethres_clin.txt
#

# Clinical + mutational signature 
python bbcar/scanmap/src/bbcar_run_scanmap.py \
    --genomic /projects/b1122/saya/scanmap_data/gene_thres_thres010_abs.pik \
    --confound /projects/b1122/saya/scanmap_data/pt_demo_clinmut_intindex.csv \
    --label /projects/b1122/saya/scanmap_data/bbcar_label_intindex.csv \
    --outdir /projects/b1122/saya/scanmap_data/out \
    --index /projects/b1122/saya/scanmap_data \
    > bbcar/scanmap/results/results_genethres_clin_mut.txt
#

# Clinical + mutational signature + somatic mutation in driver genes 
python bbcar/scanmap/src/bbcar_run_scanmap.py \
    --genomic /projects/b1122/saya/scanmap_data/gene_thres_thres010_abs.pik \
    --confound /projects/b1122/saya/scanmap_data/pt_demo_clin_mut_driversomatic_intindex.csv \
    --label /projects/b1122/saya/scanmap_data/bbcar_label_intindex.csv \
    --outdir /projects/b1122/saya/scanmap_data/out \
    --index /projects/b1122/saya/scanmap_data \
    > bbcar/scanmap/results/results_genethres_clin_mut_driversomatic.txt
#

# Clinical + PRS SNPs (removed ones appearing in <=3 patients) 
python bbcar/scanmap/src/bbcar_run_scanmap.py \
    --genomic /projects/b1122/saya/scanmap_data/gene_thres_thres010_abs.pik \
    --confound /projects/b1122/saya/scanmap_data/clin_prssnp_rm3.csv \
    --label /projects/b1122/saya/scanmap_data/bbcar_label_intindex.csv \
    --outdir /projects/b1122/saya/scanmap_data/out \
    --index /projects/b1122/saya/scanmap_data \
    > bbcar/scanmap/results/results_genethres_clin_prs.txt
#

# Clinical + mutational signature + PRS SNPs 
python bbcar/scanmap/src/bbcar_run_scanmap.py \
    --genomic /projects/b1122/saya/scanmap_data/gene_thres_thres010_abs.pik \
    --confound /projects/b1122/saya/scanmap_data/clinmut_prssnp_rm3.csv \
    --label /projects/b1122/saya/scanmap_data/bbcar_label_intindex.csv \
    --outdir /projects/b1122/saya/scanmap_data/out \
    --index /projects/b1122/saya/scanmap_data \
    > bbcar/scanmap/results/results_genethres_clin_mut_prs.txt
#

