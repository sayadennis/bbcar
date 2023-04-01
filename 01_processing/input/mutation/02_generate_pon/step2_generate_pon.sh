#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicslong
#SBATCH -n 1
#SBATCH -t 240:00:00
#SBATCH --mem=240G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="pon"
#SBATCH --output=/projects/b1131/saya/bbcar/out/generate_pon_step2.out

cd /projects/b1131/saya/bbcar/

din=/projects/b1131/saya/bbcar/data/01_alignment/germline/aligned
dout=/projects/b1131/saya/bbcar/data/02a_mutation/02_variant_calls/germline_only

interval='/projects/b1122/gannon/bbcar/RAW_data/int_lst/SureSelect_v6/hg38.preprocessed.interval_list'
ref='/projects/p30791/hg38_ref/hg38.fa'

module purge all
export PATH="/projects/b1131/saya/bbcar/tools/gatk-4.2.5.0:$PATH"
# module load gatk/4.1.0

# Instructions taken from https://gatk.broadinstitute.org/hc/en-us/articles/360035531132

# # # Step 1: Run Mutect2 in tumor-only mode for each normal sample:
# for sampleid in $(cat /projects/b1131/saya/bbcar/sample_names_germline.txt); do 
#     gatk Mutect2 -R $ref -I $din/${sampleid}_bqsr.bam --max-mnp-distance 0 -O $dout/${sampleid}_bqsr.vcf.gz;
# done
# #

# Step 2: Create a GenomicsDB from the normal Mutect2 calls:
rm -r /projects/b1042/lyglab/saya/bbcar/pon_tmp/*
rm -rf /projects/b1042/lyglab/saya/bbcar/pon_db
gatk GenomicsDBImport -R $ref -L $interval --tmp-dir /projects/b1042/lyglab/saya/bbcar/pon_tmp/ --genomicsdb-workspace-path /projects/b1042/lyglab/saya/bbcar/pon_db $(for x in $(ls -1 $dout/*_bqsr.vcf.gz); do echo -n "-V $x "; done)
#

# Step 3: Combine the normal calls using CreateSomaticPanelOfNormals:
cd /projects/b1042/lyglab/saya/bbcar/
gatk CreateSomaticPanelOfNormals \
    -V gendb://pon_db \
    -R $ref \
    --output $dout/bbcar_pon.vcf.gz
#

rm -rf /projects/b1042/lyglab/saya/bbcar/pon_db

# --germline-resource /projects/b1131/saya/bbcar/genome_resources/GATK/af-only-gnomad.hg38.vcf.gz 
