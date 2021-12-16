#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 24:00:00
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="gistic90"
#SBATCH --output=bbcar/out/03_run_gistic2_conf90.out

din=/projects/b1122/saya/02_gistic2_input
dout=/projects/b1122/saya/03_gistic2_out_conf90
refgenefile=/home/srd6051/GISTIC2/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat
# intervals="/projects/b1122/gannon/bbcar/RAW_data/int_lst/SureSelect_v6/hg38.preprocessed.interval_list"

cd /home/srd6051/GISTIC2

###############################
#### All samplies combined ####
###############################

## All - cases and controls combined 
./gistic2 -b $dout/all -seg $din/gistic_input_all.tsv -refgene $refgenefile -genegistic 0 -smallmem 0 -brlen 0.98 -conf 0.90 -armpeel 0 -savegene 1 -gcm extreme

###############################################
#### All cases and all controls separately ####
###############################################

## Cases
./gistic2 -b $dout/all_cases -seg $din/gistic_input_all_cases.tsv -refgene $refgenefile -genegistic 0 -smallmem 0 -brlen 0.98 -conf 0.90 -armpeel 0 -savegene 1 -gcm extreme
# -mk $intervals 

## Controls
./gistic2 -b $dout/all_controls -seg $din/gistic_input_all_controls.tsv -refgene $refgenefile -genegistic 0 -smallmem 0 -brlen 0.98 -conf 0.90 -armpeel 0 -savegene 1 -gcm extreme

#########################################################
#### Training cases and training controls separately ####
#########################################################

## Training cases
./gistic2 -b $dout/train_cases -seg $din/gistic_input_train_cases.tsv -refgene $refgenefile -genegistic 0 -smallmem 0 -brlen 0.98 -conf 0.90 -armpeel 0 -savegene 1 -gcm extreme

## Training controls 
./gistic2 -b $dout/train_controls -seg $din/gistic_input_train_controls.tsv -refgene $refgenefile -genegistic 0 -smallmem 0 -brlen 0.98 -conf 0.90 -armpeel 0 -savegene 1 -gcm extreme

###################################################################################
#### (Training cases + test all) and (training controls + test all) separately ####
###################################################################################

## Training cases + test all 
./gistic2 -b $dout/train_cases_test_all -seg $din/gistic_input_train_cases_test_all.tsv -refgene $refgenefile -genegistic 0 -smallmem 0 -brlen 0.98 -conf 0.90 -armpeel 0 -savegene 1 -gcm extreme

## Training controls + test all 
./gistic2 -b $dout/train_controls_test_all -seg $din/gistic_input_train_controls_test_all.tsv -refgene $refgenefile -genegistic 0 -smallmem 0 -brlen 0.98 -conf 0.90 -armpeel 0 -savegene 1 -gcm extreme
