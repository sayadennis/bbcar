#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 18:00:00
#SBATCH -N 1
#SBATCH -n 5
#SBATCH --array=0-53
#SBATCH --mem=50G
#SBATCH --job-name=bwa_ger%a
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/b1131/saya/bbcar/out/align_bwa_germline%a.out

module purge all
module load bwa
module load picard/2.6.0
module load gatk/4.1.0
module load samtools

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/bbcar/data/sample_ids_all_germline.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

set -uex

## Reference genome
FA='/projects/p30791/hg38_ref/hg38.fa'

# index the reference genome
# bwa index -p /projects/p30791/hg38_ref/hg38.fa $FA

## Software and interval locations
pic='java -jar /software/picard/2.6.0/picard-tools-2.6.0/picard.jar'
interval='/projects/b1122/gannon/bbcar/RAW_data/int_lst/SureSelect_v6/hg38.preprocessed.interval_list' # germline = V6 since all germline were seq'ed at Indiana
# interval='/projects/b1122/gannon/bbcar/RAW_data/int_lst/SureSelect_v5/hg38/hg38.preprocessed.interval_list' # V5
dbsnp='/projects/b1131/saya/bbcar/genome_resources/Homo_sapiens_assembly38.dbsnp138.vcf' # is this file appropriate for this?
gold1000Indel='/projects/b1131/saya/bbcar/genome_resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz' # is this file appropriate for this?

## Output directories
metrics='/projects/b1131/saya/bbcar/data/01_alignment/germline/metrics'
rec_tables='/projects/b1131/saya/bbcar/data/01_alignment/germline/recal_tables'
interim='/projects/b1131/saya/bbcar/data/01_alignment/germline/interim'
aligned='/projects/b1131/saya/bbcar/data/01_alignment/germline/aligned'

cd /projects/b1131/saya/bbcar/data/00_raw/germline/${sampleid}/

for fpattern in $(awk -F'[_\.]' '{print $1"_"$2"_"$3"_"$4"_"$5"_"$6}' <(ls *.fastq.gz) 2>/dev/null | cut -d '_' -f 1-3 | sort -u); do 
    # echo $fpattern

    R1=/projects/b1131/saya/bbcar/data/00_raw/germline/${sampleid}/${fpattern}_R1.fastq.gz
    R2=/projects/b1131/saya/bbcar/data/00_raw/germline/${sampleid}/${fpattern}_R2.fastq.gz

    header=$(zcat $R1 | head -n 1)
    id=$(echo $header | head -n 1 | cut -f 3-4 -d':' | sed 's/@//' | sed 's/:/_/g')
    echo "@RG\tID:${id}\tSM:${sampleid}_germline\tLB:library1\tPL:ILLUMINA\tPU:${id}"

    bwa mem -M -t 24 -R $(echo "@RG\tID:${id}\tSM:${sampleid}_germline\tLB:library1\tPL:ILLUMINA\tPU:${id}") $FA $R1 $R2 > $interim/${fpattern}.sam

    $pic SortSam I=${interim}/${fpattern}.sam O=${interim}/${fpattern}_sorted.bam SORT_ORDER=coordinate
done

$pic MergeSamFiles CREATE_INDEX=true $(for x in $(ls -1 ${interim}/${sampleid}_*_sorted.bam); do echo -n "I=$x "; done) O=${interim}/${sampleid}_final_sorted.bam USE_THREADING=true

rm -f ${interim}/${sampleid}_*_*_.sam
rm -f ${interim}/${sampleid}_*_*_*_sorted.bam

$pic MarkDuplicates I=${interim}/${sampleid}_final_sorted.bam O=${interim}/${sampleid}_dup.bam M=${metrics}/${sampleid}_germline_reads.mdup.metrics.txt

rm -f ${interim}/${sampleid}_final_sorted.bam
rm -f ${interim}/${sampleid}_final_sorted.bai

gatk BaseRecalibrator -I ${interim}/${sampleid}_dup.bam -R $FA --known-sites $dbsnp --known-sites $gold1000Indel -O ${rec_tables}/${sampleid}_germline_recal_data.table

samtools index ${interim}/${sampleid}_dup.bam

gatk ApplyBQSR -R $FA -I ${interim}/${sampleid}_dup.bam --bqsr-recal-file ${rec_tables}/${sampleid}_germline_recal_data.table -L ${interval} -O ${aligned}/${sampleid}_bqsr.bam

rm -f ${interim}/${sampleid}_dup.bam
rm -f ${interim}/${sampleid}_dup.bam.bai
