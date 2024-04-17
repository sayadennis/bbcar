#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics-himem
#SBATCH -n 1
#SBATCH -t 168:00:00
#SBATCH --mem=500G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="pon_beta"
#SBATCH --output=bbcar/out/generate_pon_bbcar_beta.out

cd /projects/b1131/saya/bbcar/

din=/projects/b1131/saya/bbcar/01_alignment/germline/aligned
dout=/projects/b1131/saya/bbcar/02_variant_calls/germline_only

interval='/projects/b1122/gannon/bbcar/RAW_data/int_lst/SureSelect_v6/hg38.preprocessed.interval_list'
ref='/projects/p30791/hg38_ref/hg38.fa'

module purge all
export PATH="/projects/b1131/saya/bbcar/tools/gatk-4.2.5.0:$PATH"
# module load gatk/4.1.0

# Instructions taken from https://gatk.broadinstitute.org/hc/en-us/articles/360035531132

# ## Step 1: Run Mutect2 in tumor-only mode for each normal sample:
# for sampleid in $(cat /projects/b1131/saya/bbcar/sample_names_germline.txt); do 
#     gatk Mutect2 -R $ref -I $din/${sampleid}_bqsr.bam --max-mnp-distance 0 -O $dout/${sampleid}_bqsr.vcf.gz;
# done
# #

## Step 2: Create an argument file 
touch /projects/b1131/saya/bbcar/normals_for_pon_vcf.args # create file if it doesn't exist 
> /projects/b1131/saya/bbcar/normals_for_pon_vcf.args # empty contents of the file 
ls ${dout}/*_bqsr.vcf.gz > /projects/b1131/saya/bbcar/normals_for_pon_vcf.args

## Step 3: Combine the normal calls using CreateSomaticPanelOfNormals:
gatk CreateSomaticPanelOfNormals \
   -vcfs /projects/b1131/saya/bbcar/normals_for_pon_vcf.args \
   -O /projects/b1131/saya/bbcar/bbcar_pon.vcf.gz
#

# --germline-resource /projects/b1131/saya/bbcar/genome_resources/GATK/af-only-gnomad.hg38.vcf.gz 
