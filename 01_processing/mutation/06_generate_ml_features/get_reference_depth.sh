#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 48:00:00
#SBATCH --mem=10G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=refdp
#SBATCH --output=bbcar/out/get_reference_depth.out

module purge all
module load bwa
module load picard/2.6.0
module load gatk/4.1.0
module load samtools/1.10.1

ref='/projects/p30791/hg38_ref/hg38.fa'

###################
#### Alignment ####
###################

cd /projects/p30791/

set -uex

pic='java -jar /software/picard/2.6.0/picard-tools-2.6.0/picard.jar'
metrics=/projects/b1122/saya/bbcar/dup_metrics
# rec_tables=/projects/b1122/saya/bbcar/recal_tables
# interval='/projects/b1122/Zexian/tools/DNAtools/S07604514_Padded.bed'

bwa mem -M -t 24 -R $ref $ref > /projects/p30791/hg38_ref/hg38.sam # this can't be right lol

$pic SortSam I=/projects/p30791/hg38_ref/hg38.sam O=/projects/p30791/hg38_ref/hg38.bam SORT_ORDER=coordinate

rm -f /projects/p30791/hg38_ref/hg38.sam
# rm -f hg38.bam

$pic MarkDuplicates I=/projects/p30791/hg38_ref/hg38.bam O=/projects/p30791/hg38_ref/hg38_dup.bam M=$metrics/hg38_reads.mdup.metrics.txt

# gatk BaseRecalibrator -I 1419_dup.bam -R $FA --known-sites $dbsnp --known-sites $gold1000Indel -O $rec_tables/1419_tissue_recal_data.table

samtools index hg38_dup.bam

############################
#### Calculate coverage ####
############################

samtools depth /projects/p30791/hg38_ref/deduped_MA605.bam > /projects/p30791/hg38_ref/deduped_MA605.coverage

################################################
#### Select depth of the region of interest ####
################################################

## Step 2

# To select the coverage for a particular chromosome (Chr#1 in my case)
awk '$1 == 1 {print $0}' deduped_MA605.coverage > chr1_MA605.coverage

# # If the chrosomosome has string characters it can be adjusted as
# awk '$1 == "chr2" {print $0}' deduped_MA605.coverage > chr2_MA605.coverage

# # Step #2) can be done in one step like this
# gawk '{/^[0-9]/{print >$1".coverag"}' ./deduped_MA605.coverag

