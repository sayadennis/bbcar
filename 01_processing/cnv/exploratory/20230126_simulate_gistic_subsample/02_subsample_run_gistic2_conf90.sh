#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --array=0-4
#SBATCH --mem=1G
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="subgistic%a"
#SBATCH --output=/projects/b1131/saya/bbcar/out/subsample_run_gistic2_conf90%a.out

din="/projects/b1131/saya/bbcar/data/02b_cnv/20230126_simulate_gistic_subsample/gistic_input"
dout="/projects/b1131/saya/bbcar/data/02b_cnv/20230126_simulate_gistic_subsample/gistic_output"
refgenefile="/projects/p30791/GISTIC2/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat"
# intervals="/projects/b1122/gannon/bbcar/RAW_data/int_lst/SureSelect_v6/hg38.preprocessed.interval_list"

cd /projects/p30791/GISTIC2/

mkdir -p $dout/subset${SLURM_ARRAY_TASK_ID}

./gistic2 -b $dout/subset${SLURM_ARRAY_TASK_ID} -seg $din/gistic_input_subset${SLURM_ARRAY_TASK_ID}.tsv -refgene $refgenefile -genegistic 0 -smallmem 0 -brlen 0.98 -conf 0.90 -armpeel 0 -savegene 1 -gcm extreme
