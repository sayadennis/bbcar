#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --mem=50G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/home/srd6051/bbcar/out/bwa_1419.out
#SBATCH -J bwa_1419

module purge all
module load bwa
module load picard/2.6.0
module load gatk/4.1.0
module load samtools

cd /projects/p30791/bbcar_test

set -uex

FA=/projects/b1122/gannon/ref_seqs/ucsc_hg38/hg38.fa
# Zex_germ=/projects/b1122/Zexian/Alignment/Germline_37/RAW_data
# New_germ=/projects/b1122/Zexian/Alignment/BBCAR_GC
# tissue_fq=/projects/b1122/Zexian/Alignment/BBCAR/RAW_data
# idx=/projects/p30007/gannon/bbcar/bwa_algn/refs
pic='java -jar /software/picard/2.6.0/picard-tools-2.6.0/picard.jar'
metrics=/projects/b1122/saya/bbcar/dup_metrics
# gold1000Indel=/projects/b1122/Zexian/reference/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
## ^FIGURE THIS OUT --DO WE HAVE HG38 VERSION AND/OR IS THIS SPECIFICATION OF KNOWN SITES NECESSARY?
# dbsnp=/projects/b1122/Zexian/reference/hg19/dbsnp_138.hg19.vcf
## ^same with this. Find a hg38 version.
rec_tables=/projects/b1122/saya/bbcar/recal_tables
interval='/projects/b1122/Zexian/tools/DNAtools/S07604514_Padded.bed'

R1=/projects/b1122/saya/raw/bbb_tissue/1419/1419_S36_L005_R2.fastq.gz
R2=/projects/b1122/saya/raw/bbb_tissue/1419/1419_S36_L005_R1.fastq.gz

header=$(zcat $R1 | head -n 1)
id=$(echo $header | head -n 1 | cut -f 3-4 -d":" | sed 's/@//' | sed 's/:/_/g')
echo "@RG\tID:$id\tSM:1419_germ\tLB:library1\tPL:ILLUMINA\tPU:$id"

bwa mem -M -t 24 -R $(echo "@RG\tID:$id\tSM:1419_germ\tLB:library1\tPL:ILLUMINA\tPU:$id") $FA $R1 $R2 > 1419_1.sam

$pic SortSam I=1419_1.sam O=1419_1_sorted.bam SORT_ORDER=coordinate

R1=/projects/b1122/saya/raw/bbb_tissue/1419/1419_S36_L008_R2.fastq.gz
R2=/projects/b1122/saya/raw/bbb_tissue/1419/1419_S36_L008_R1.fastq.gz

header=$(zcat $R1 | head -n 1)
id=$(echo $header | head -n 1 | cut -f 3-4 -d":" | sed 's/@//' | sed 's/:/_/g')
echo "@RG\tID:$id\tSM:1419_germ\tLB:library1\tPL:ILLUMINA\tPU:$id"

bwa mem -M -t 24 -R $(echo "@RG\tID:$id\tSM:1419_germ\tLB:library1\tPL:ILLUMINA\tPU:$id") $FA $R1 $R2 > 1419_2.sam

$pic SortSam I=1419_2.sam O=1419_2_sorted.bam SORT_ORDER=coordinate

R1=/projects/b1122/saya/raw/bbb_tissue/1419/1419_S36_L007_R1.fastq.gz
R2=/projects/b1122/saya/raw/bbb_tissue/1419/1419_S36_L007_R2.fastq.gz

header=$(zcat $R1 | head -n 1)
id=$(echo $header | head -n 1 | cut -f 3-4 -d":" | sed 's/@//' | sed 's/:/_/g')
echo "@RG\tID:$id\tSM:1419_germ\tLB:library1\tPL:ILLUMINA\tPU:$id"

bwa mem -M -t 24 -R $(echo "@RG\tID:$id\tSM:1419_germ\tLB:library1\tPL:ILLUMINA\tPU:$id") $FA $R1 $R2 > 1419_3.sam

$pic SortSam I=1419_3.sam O=1419_3_sorted.bam SORT_ORDER=coordinate

R1=/projects/b1122/saya/raw/bbb_tissue/1419/1419_S36_L006_R2.fastq.gz
R2=/projects/b1122/saya/raw/bbb_tissue/1419/1419_S36_L006_R1.fastq.gz

header=$(zcat $R1 | head -n 1)
id=$(echo $header | head -n 1 | cut -f 3-4 -d":" | sed 's/@//' | sed 's/:/_/g')
echo "@RG\tID:$id\tSM:1419_germ\tLB:library1\tPL:ILLUMINA\tPU:$id"

bwa mem -M -t 24 -R $(echo "@RG\tID:$id\tSM:1419_germ\tLB:library1\tPL:ILLUMINA\tPU:$id") $FA $R1 $R2 > 1419_4.sam

$pic SortSam I=1419_4.sam O=1419_4_sorted.bam SORT_ORDER=coordinate


$pic MergeSamFiles CREATE_INDEX=true I=1419_1_sorted.bam I=1419_2_sorted.bam I=1419_3_sorted.bam I=1419_4_sorted.bam O=1419_final_sorted.bam USE_THREADING=true

rm -f 1419_1.sam
rm -f 1419_2.sam
rm -f 1419_3.sam
rm -f 1419_4.sam
rm -f 1419_1_sorted.bam
rm -f 1419_2_sorted.bam
rm -f 1419_3_sorted.bam
rm -f 1419_4_sorted.bam

$pic MarkDuplicates I=1419_final_sorted.bam O=1419_dup.bam M=$metrics/1419_tissue_reads.mdup.metrics.txt

rm -f 1419_final_sorted.bam
rm -f 1419_final_sorted.bai

gatk BaseRecalibrator -I 1419_dup.bam -R $FA --known-sites $dbsnp --known-sites $gold1000Indel -O $rec_tables/1419_tissue_recal_data.table

samtools index 1419_dup.bam

gatk ApplyBQSR -R $FA -I 1419_dup.bam --bqsr-recal-file $rec_tables/1419_tissue_recal_data.table -L $interval -O 1419t_bqsr.bam

rm -f 1419_dup.bam
rm -f 1419_dup.bam.bai
