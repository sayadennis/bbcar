#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --array=0-53
#SBATCH -n 1
#SBATCH -t 4:00:00
#SBATCH --mem=3G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="haplo_\${SLURM_ARRAY_TASK_ID}"
#SBATCH --output=bbcar/out/call_variants_germline_haplotype_%a.out

cd /projects/b1131/saya/bbcar/

## Load GATK 
module purge all
export PATH="/projects/b1131/saya/bbcar/tools/gatk-4.2.5.0:$PATH"

#### Create sample names TXT file for job array (RUN THIS IN ADVANCE) ####
# ls /projects/b1131/saya/bbcar/01_alignment/tissue/aligned/*.bam | xargs -n1 basename | tr "_" "\n" | grep "t" > /projects/b1131/saya/bbcar/sample_names_tumor_only.txt

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/bbcar/sample_names_germline.txt

## Set input and output directories
din='/projects/b1131/saya/bbcar/01_alignment/germline/aligned'
dout='/projects/b1131/saya/bbcar/02_variant_calls/germline_only'

## Set reference, interval, germline resource, and PON filenames 
ref='/projects/p30791/hg38_ref/hg38.fa'
interval='/projects/b1122/gannon/bbcar/RAW_data/int_lst/SureSelect_v6/hg38.preprocessed.interval_list'
germres='/projects/b1131/saya/bbcar/genome_resources/GATK/af-only-gnomad.hg38.vcf.gz'

## Instructions taken from https://gatk.broadinstitute.org/hc/en-us/articles/360035531132 
## Under section "A step-by-step guide to the new Mutect2 Read Orientation Artifacts Workflow"

## Create output with raw data used to learn the orientation bias model
gatk --java-options "-Xmx3g" HaplotypeCaller \
        -R $ref \
        --intervals $interval \
        -I $din/${input_args[$SLURM_ARRAY_TASK_ID]}_bqsr.bam \
        -O $dout/${input_args[$SLURM_ARRAY_TASK_ID]}_haplotype.vcf
#
