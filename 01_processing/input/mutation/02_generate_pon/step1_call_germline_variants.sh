#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -n 1
#SBATCH --array=0-53
#SBATCH -t 3:00:00
#SBATCH --mem=3G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --job-name="varcallpon"
#SBATCH --output=/projects/b1131/saya/bbcar/out/generate_pon_step1_call_variants.out

cd /projects/b1131/saya/bbcar/

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/bbcar/data/sample_ids_all_germline.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

din=/projects/b1131/saya/bbcar/data/01_alignment/germline/aligned
dout=/projects/b1131/saya/bbcar/data/02a_mutation/02_variant_calls/germline_only

interval='/projects/b1122/gannon/bbcar/RAW_data/int_lst/SureSelect_v6/hg38.preprocessed.interval_list'
ref='/projects/p30791/hg38_ref/hg38.fa'

module purge all
export PATH="/projects/b1131/saya/bbcar/tools/gatk-4.2.5.0:$PATH"
# module load gatk/4.1.0

# Instructions taken from https://gatk.broadinstitute.org/hc/en-us/articles/360035531132

# # Step 1: Run Mutect2 in tumor-only mode for each normal sample:
gatk Mutect2 -R $ref -I ${din}/${sampleid}_bqsr.bam --max-mnp-distance 0 -O ${dout}/${sampleid}_bqsr.vcf.gz
#
