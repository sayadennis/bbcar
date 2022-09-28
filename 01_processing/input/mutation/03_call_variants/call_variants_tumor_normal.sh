#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --array=0-51
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --mem=5G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="mutect2_tn_%a"
#SBATCH --output=bbcar/out/call_variants_tumor_normal_%a.out

cd /projects/b1131/saya/bbcar/

## Load GATK 
module purge all
export PATH="/projects/b1131/saya/bbcar/tools/gatk-4.2.5.0:$PATH"

################################
## Create sample names TXT file for job array -- run this in Python for ease of handling ##
# import os
# germ = os.listdir('/projects/b1131/saya/bbcar/01_alignment/germline/aligned')
# germ = [x.split('_')[0] for x in germ]
# tiss = os.listdir('/projects/b1131/saya/bbcar/01_alignment/tissue/aligned')
# tiss = [x.split('t')[0] for x in tiss]
# both = []
# for item in germ:
#       if item in tiss:
#               if item not in both:
#                       both.append(item)
# 
# with open('/projects/b1131/saya/bbcar/sample_names_tumor_normal.txt', 'w') as f:
#       for item in both:
#               f.write(f'{item}\n')
################################

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < sample_names_tumor_normal.txt

## Set input and output directories
din_tiss='/projects/b1131/saya/bbcar/01_alignment/tissue/aligned'
din_germ='/projects/b1131/saya/bbcar/01_alignment/germline/aligned'
dout='/projects/b1131/saya/bbcar/02_variant_calls/tumor_normal'

## Set reference, interval, germline resource, and PON filenames 
ref='/projects/p30791/hg38_ref/hg38.fa'
interval='/projects/b1122/gannon/bbcar/RAW_data/int_lst/SureSelect_v6/hg38.preprocessed.interval_list'
germres='/projects/b1131/saya/bbcar/genome_resources/GATK/af-only-gnomad.hg38.vcf.gz'

## Instructions taken from https://gatk.broadinstitute.org/hc/en-us/articles/360035531132 
## Under section "A step-by-step guide to the new Mutect2 Read Orientation Artifacts Workflow"

# sampleid='1004t'

##########################
#### With generic PON ####
##########################

pon='/projects/b1131/saya/bbcar/genome_resources/GATK/1000g_pon.hg38.vcf.gz'

## Create output with raw data used to learn the orientation bias model
gatk Mutect2 \
        -R $ref \
        -L $interval \
        -I $din_tiss/${input_args[$SLURM_ARRAY_TASK_ID]}t_bqsr.bam \
        -I $din_germ/${input_args[$SLURM_ARRAY_TASK_ID]}_bqsr.bam \
        --normal-sample ${input_args[$SLURM_ARRAY_TASK_ID]}_germline \
        -germline-resource $germres \
        -pon $pon \
        --f1r2-tar-gz $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_f1r2.tar.gz \
        -O $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_unfiltered.vcf
#

## Pass this raw data to LearnReadOrientationModel
gatk LearnReadOrientationModel -I $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_f1r2.tar.gz -O $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_read-orientation-model.tar.gz

## Run GetPileupSummaries to summarize read support for a set number of known variant sites.
gatk GetPileupSummaries \
    --input $din/${input_args[$SLURM_ARRAY_TASK_ID]}_bqsr.bam \
    --variant $germres \
    --intervals $interval \
    --output $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_getpileupsummaries.table
#

## Estimate contamination with CalculateContamination.
gatk CalculateContamination \
        -I $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_getpileupsummaries.table \
        -tumor-segmentation $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_segments.table \
        -O $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_calculatecontamination.table
#

## Finally, pass the learned read orientation model to FilterMutectCallswith the -ob-priors argument:
gatk FilterMutectCalls -V $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_unfiltered.vcf \
        -R $ref \
        --ob-priors $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_read-orientation-model.tar.gz \
        -O $dout/${input_args[$SLURM_ARRAY_TASK_ID]}_filtered.vcf
        # [--tumor-segmentation $dout/interim/segments.table] \
        # [--contamination-table contamination.table] \
#

########################
#### With BBCAR PON ####
########################

pon='/projects/b1131/saya/bbcar/02_variant_calls/germline_only/bbcar_pon.vcf.gz'

## Create output with raw data used to learn the orientation bias model
gatk Mutect2 \
        -R $ref \
        -L $interval \
        -I $din_tiss/${input_args[$SLURM_ARRAY_TASK_ID]}t_bqsr.bam \
        -I $din_germ/${input_args[$SLURM_ARRAY_TASK_ID]}_bqsr.bam \
        --normal-sample ${input_args[$SLURM_ARRAY_TASK_ID]}_germline \
        -germline-resource $germres \
        -pon $pon \
        --f1r2-tar-gz $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_f1r2_bbcarpon.tar.gz \
        -O $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_unfiltered_bbcarpon.vcf
#

## Pass this raw data to LearnReadOrientationModel
gatk LearnReadOrientationModel -I $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_f1r2_bbcarpon.tar.gz -O $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_read-orientation-model_bbcarpon.tar.gz

## Run GetPileupSummaries to summarize read support for a set number of known variant sites.
gatk GetPileupSummaries \
    --input $din/${input_args[$SLURM_ARRAY_TASK_ID]}_bqsr.bam \
    --variant $germres \
    --intervals $interval \
    --output $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_getpileupsummaries_bbcarpon.table
#

## Estimate contamination with CalculateContamination.
gatk CalculateContamination \
        -I $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_getpileupsummaries_bbcarpon.table \
        -tumor-segmentation $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_segments_bbcarpon.table \
        -O $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_calculatecontamination_bbcarpon.table
#

## Finally, pass the learned read orientation model to FilterMutectCallswith the -ob-priors argument:
gatk FilterMutectCalls -V $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_unfiltered_bbcarpon.vcf \
        -R $ref \
        --ob-priors $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_read-orientation-model_bbcarpon.tar.gz \
        -O $dout/${input_args[$SLURM_ARRAY_TASK_ID]}_filtered_bbcarpon.vcf
        # [--tumor-segmentation $dout/interim/segments.table] \
        # [--contamination-table contamination.table] \
#

