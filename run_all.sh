#!/bin/bash

cd $HOME

# Source configuration file
source bbcar/repo/run_all.config.ini

# Export all config variables so they can be accessed by the submitted jobs
export PROJECT_DIR
export RAW_DATA_DIR
export GENOME_RESOURCES_DIR
export REF_GENOME_FA
export INTLIST_V5
export INTLIST_V6
export TISSUE_SAMPLEID_FILE
export GERMLINE_SAMPLEID_FILE
export V5_SAMPLEID_FILE
export LABEL_FILE

ARRAY_MAX_TISSUE=$(( $(wc -l < ${TISSUE_SAMPLEID_FILE}) - 1))
ARRAY_MAX_GERMLINE=$(( $(wc -l < ${GERMLINE_SAMPLEID_FILE}) - 1))

# Step 1: alignment on tissue samples
jid0=($(sbatch --array=0-${ARRAY_MAX_TISSUE} bbcar/repo/01_processing/input/mutation/01_align/align_tissues.sh))

echo "jid0 ${jid0[-1]}" >> slurm_ids

# Step 2: alignment on germline samples
jid1=($(sbatch --array=0-${ARRAY_MAX_GERMLINE} --dependency=afterok:${jid0[-1]} --export=DEPENDENTJOB=${jid0[-1]} bbcar/repo/01_processing/input/mutation/01_align/align_germlines.sh))

echo "jid1 ${jid1[-1]}" >> slurm_ids

#jid2=($(sbatch --dependency=afterok:${jid1[-1]} --export=DEPENDENTJOB=${jid1[-1]} example_submit.sh))

#echo "jid2 ${jid2[-1]}" >> slurm_ids
