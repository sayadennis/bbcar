#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomicslong
#SBATCH -n 1
#SBATCH -t 240:00:00
#SBATCH --mem=240G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="pon"
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/generate_pon_step2.out

cd /projects/b1042/lyglab/saya/

din='/projects/b1131/saya/new_bbcar/data/01_alignment/germline/aligned'
dout='/projects/b1131/saya/new_bbcar/data/02a_mutation/02_variant_calls/germline_only'

interval='/projects/b1122/gannon/bbcar/RAW_data/int_lst/SureSelect_v6/hg38.preprocessed.interval_list'
ref='/projects/p30791/hg38_ref/hg38.fa'

module purge all
module load singularity

gatk() {
    singularity exec -B /projects:/projects /projects/p30791/gatk_4.5.0.sif gatk "$@"
}

# Instructions taken from https://gatk.broadinstitute.org/hc/en-us/articles/360035531132

# Step 2: Create a GenomicsDB from the normal Mutect2 calls:
PON_TMP='/projects/b1042/lyglab/saya/new_bbcar/pon_tmp'
PON_DB='/projects/b1042/lyglab/saya/new_bbcar/pon_db'

mkdir -p ${PON_TMP}
mkdir -p ${PON_DB}

rm -r /projects/b1042/lyglab/saya/new_bbcar/pon_tmp/*
rm -rf /projects/b1042/lyglab/saya/new_bbcar/pon_db

gatk GenomicsDBImport \
    -R $ref \
    -L $interval \
    --tmp-dir ${PON_TMP}/ \
    --genomicsdb-workspace-path ${PON_DB} \
    $(for x in $(ls -1 $dout/*_bqsr.vcf.gz); do echo -n "-V $x "; done)

# Step 3: Combine the normal calls using CreateSomaticPanelOfNormals:
cd /projects/b1042/lyglab/saya/new_bbcar/
gatk CreateSomaticPanelOfNormals \
    -V gendb://pon_db \
    -R $ref \
    --output $dout/bbcar_pon.vcf.gz

rm -rf /projects/b1042/lyglab/saya/bbcar/pon_db

# --germline-resource /projects/b1131/saya/bbcar/genome_resources/GATK/af-only-gnomad.hg38.vcf.gz 
