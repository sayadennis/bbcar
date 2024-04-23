#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -n 1
#SBATCH --array=52,53,54,57,58,59,60,61,62,63,64,65,66,67,69,70,71,72,73,74,77,78,79,80,81,82,83,84,85,86,87,88,89,90,92
#SBATCH -t 12:00:00
#SBATCH --mem=3G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --job-name="varcallpon"
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/generate_pon_step1_call_variants.out

cd /projects/b1131/saya/new_bbcar/

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_germline.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

din=/projects/b1131/saya/new_bbcar/data/01_alignment/germline/aligned
dout=/projects/b1131/saya/new_bbcar/data/02a_mutation/02_variant_calls/germline_only

mkdir -p $dout

interval='/projects/b1122/gannon/bbcar/RAW_data/int_lst/SureSelect_v6/hg38.preprocessed.interval_list'
ref='/projects/p30791/hg38_ref/hg38.fa'

module purge all
module load singularity

gatk() {
    singularity exec -B /projects:/projects /projects/p30791/gatk_4.5.0.sif gatk "$@"
}

# Instructions taken from https://gatk.broadinstitute.org/hc/en-us/articles/360035531132

# Step 1: Run Mutect2 in tumor-only mode for each normal sample:
gatk Mutect2 -R $ref -I ${din}/${sampleid}_bqsr.bam --max-mnp-distance 0 -O ${dout}/${sampleid}_bqsr.vcf.gz

