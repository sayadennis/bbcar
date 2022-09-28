#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --mem=50G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/home/srd6051/bbcar/out/bwa_743_tissue.out
#SBATCH --job-name=bwa_743_tissue

module purge all
module load bwa
module load picard/2.6.0
module load gatk/4.1.0
module load samtools

cd /projects/p30791/bbcar_test

set -uex

## Reference genome
FA='/projects/p30791/hg38_ref/hg38.fa'

# index the reference genome
# bwa index -p /projects/p30791/hg38_ref/hg38.fa $FA

## Software and interval locations
pic='java -jar /software/picard/2.6.0/picard-tools-2.6.0/picard.jar'
interval='/projects/b1122/Zexian/tools/DNAtools/S07604514_Padded.bed'
dbsnp='/projects/b1131/saya/bbcar/dbsnp/Homo_sapiens_assembly38.dbsnp138.vcf' # is this file appropriate for this?
gold1000Indel='/projects/b1131/saya/bbcar/dbsnp/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz' # is this file appropriate for this?

## Output directories
metrics='/projects/b1131/saya/bbcar/01_alignment/tissue/metrics'
rec_tables='/projects/b1131/saya/bbcar/01_alignment/tissue/recal_tables'
interim='/projects/b1131/saya/bbcar/01_alignment/tissue/interim'
aligned='/projects/b1131/saya/bbcar/01_alignment/tissue/aligned'

R1=/projects/b1131/saya/bbcar/raw/tissue/743/743_S129_L008_R1.fastq.gz
R2=/projects/b1131/saya/bbcar/raw/tissue/743/743_S129_L008_R2.fastq.gz

header=$(zcat $R1 | head -n 1)
id=$(echo $header | head -n 1 | cut -f 3-4 -d':' | sed 's/@//' | sed 's/:/_/g')
echo '@RG\tID:$id\tSM:743_tissue\tLB:library1\tPL:ILLUMINA\tPU:$id'

bwa mem -M -t 24 -R $(echo '@RG\tID:$id\tSM:743_tissue\tLB:library1\tPL:ILLUMINA\tPU:$id') $FA $R1 $R2 > $interim/743_1.sam

$pic SortSam I=$interim/743_1.sam O=$interim/743_1_sorted.bam SORT_ORDER=coordinate

R1=/projects/b1131/saya/bbcar/raw/tissue/743/743_S129_L008_R1.fastq.gz
R2=/projects/b1131/saya/bbcar/raw/tissue/743/743_R1.fastq.gz

header=$(zcat $R1 | head -n 1)
id=$(echo $header | head -n 1 | cut -f 3-4 -d':' | sed 's/@//' | sed 's/:/_/g')
echo '@RG\tID:$id\tSM:743_tissue\tLB:library1\tPL:ILLUMINA\tPU:$id'

bwa mem -M -t 24 -R $(echo '@RG\tID:$id\tSM:743_tissue\tLB:library1\tPL:ILLUMINA\tPU:$id') $FA $R1 $R2 > $interim/743_2.sam

$pic SortSam I=$interim/743_2.sam O=$interim/743_2_sorted.bam SORT_ORDER=coordinate

R1=/projects/b1131/saya/bbcar/raw/tissue/743/743_S129_L007_R1.fastq.gz
R2=/projects/b1131/saya/bbcar/raw/tissue/743/743_S129_L007_R2.fastq.gz

header=$(zcat $R1 | head -n 1)
id=$(echo $header | head -n 1 | cut -f 3-4 -d':' | sed 's/@//' | sed 's/:/_/g')
echo '@RG\tID:$id\tSM:743_tissue\tLB:library1\tPL:ILLUMINA\tPU:$id'

bwa mem -M -t 24 -R $(echo '@RG\tID:$id\tSM:743_tissue\tLB:library1\tPL:ILLUMINA\tPU:$id') $FA $R1 $R2 > $interim/743_3.sam

$pic SortSam I=$interim/743_3.sam O=$interim/743_3_sorted.bam SORT_ORDER=coordinate

R1=/projects/b1131/saya/bbcar/raw/tissue/743/743_S35_L005_R2.fastq.gz
R2=/projects/b1131/saya/bbcar/raw/tissue/743/743_S35_L005_R1.fastq.gz

header=$(zcat $R1 | head -n 1)
id=$(echo $header | head -n 1 | cut -f 3-4 -d':' | sed 's/@//' | sed 's/:/_/g')
echo '@RG\tID:$id\tSM:743_tissue\tLB:library1\tPL:ILLUMINA\tPU:$id'

bwa mem -M -t 24 -R $(echo '@RG\tID:$id\tSM:743_tissue\tLB:library1\tPL:ILLUMINA\tPU:$id') $FA $R1 $R2 > $interim/743_4.sam

$pic SortSam I=$interim/743_4.sam O=$interim/743_4_sorted.bam SORT_ORDER=coordinate

R1=/projects/b1131/saya/bbcar/raw/tissue/743/743_S129_L006_R2.fastq.gz
R2=/projects/b1131/saya/bbcar/raw/tissue/743/743_S129_L006_R1.fastq.gz

header=$(zcat $R1 | head -n 1)
id=$(echo $header | head -n 1 | cut -f 3-4 -d':' | sed 's/@//' | sed 's/:/_/g')
echo '@RG\tID:$id\tSM:743_tissue\tLB:library1\tPL:ILLUMINA\tPU:$id'

bwa mem -M -t 24 -R $(echo '@RG\tID:$id\tSM:743_tissue\tLB:library1\tPL:ILLUMINA\tPU:$id') $FA $R1 $R2 > $interim/743_5.sam

$pic SortSam I=$interim/743_5.sam O=$interim/743_5_sorted.bam SORT_ORDER=coordinate

R1=/projects/b1131/saya/bbcar/raw/tissue/743/743_S129_L005_R2.fastq.gz
R2=/projects/b1131/saya/bbcar/raw/tissue/743/743_S129_L005_R1.fastq.gz

header=$(zcat $R1 | head -n 1)
id=$(echo $header | head -n 1 | cut -f 3-4 -d':' | sed 's/@//' | sed 's/:/_/g')
echo '@RG\tID:$id\tSM:743_tissue\tLB:library1\tPL:ILLUMINA\tPU:$id'

bwa mem -M -t 24 -R $(echo '@RG\tID:$id\tSM:743_tissue\tLB:library1\tPL:ILLUMINA\tPU:$id') $FA $R1 $R2 > $interim/743_6.sam

$pic SortSam I=$interim/743_6.sam O=$interim/743_6_sorted.bam SORT_ORDER=coordinate


$pic MergeSamFiles CREATE_INDEX=true I=$interim/743_1_sorted.bam I=$interim/743_2_sorted.bam I=$interim/743_3_sorted.bam I=$interim/743_4_sorted.bam I=$interim/743_5_sorted.bam I=$interim/743_6_sorted.bam O=$interim/743_final_sorted.bam USE_THREADING=true

rm -f $interim/743_1.sam
rm -f $interim/743_2.sam
rm -f $interim/743_3.sam
rm -f $interim/743_4.sam
rm -f $interim/743_5.sam
rm -f $interim/743_6.sam
rm -f $interim/743_1_sorted.bam
rm -f $interim/743_2_sorted.bam
rm -f $interim/743_3_sorted.bam
rm -f $interim/743_4_sorted.bam
rm -f $interim/743_5_sorted.bam
rm -f $interim/743_6_sorted.bam

$pic MarkDuplicates I=$interim/743_final_sorted.bam O=$interim/743_dup.bam M=$metrics/743_tissue_reads.mdup.metrics.txt

rm -f $interim/743_final_sorted.bam
rm -f $interim/743_final_sorted.bai

gatk BaseRecalibrator -I $interim/743_dup.bam -R $FA --known-sites $dbsnp --known-sites $gold1000Indel -O $rec_tables/743_tissue_recal_data.table

samtools index $interim/743_dup.bam

gatk ApplyBQSR -R $FA -I $interim/743_dup.bam --bqsr-recal-file $rec_tables/743_tissue_recal_data.table -L $interval -O $aligned/743_bqsr.bam

rm -f $interim/743_dup.bam
rm -f $interim/743_dup.bam.bai
