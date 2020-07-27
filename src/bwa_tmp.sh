#!/bin/bash              
[settings]

module purge all
module load samtools
module load bwa
module load bowtie2/2.2.6
module load python
module load R
module load picard/2.6.0
module load gatk/4.1.0

cd /projects/b1042/ClareLab/GATK_PON/bam/germ

set -uex

FA=/projects/p30007/Zexian/reference/hg19/ucsc.hg19.fasta
Zex_germ=/projects/p30007/Zexian/Alignment/Germline_37/RAW_data
New_germ=/projects/p30007/Zexian/Alignment/BBCAR_GC
tissue_fq=/projects/p30007/Zexian/Alignment/BBCAR/RAW_data
idx=/projects/p30007/gannon/bbcar/bwa_algn/refs
pic='java -jar /software/picard/2.6.0/picard-tools-2.6.0/picard.jar'
metrics=/projects/p30007/gannon/bbcar/dup_metrics
gold1000Indel=/projects/p30007/Zexian/reference/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
dbsnp=/projects/p30007/Zexian/reference/hg19/dbsnp_138.hg19.vcf
rec_tables=/projects/p30007/gannon/bbcar/recal_table
interval='/projects/p30007/Zexian/tools/DNAtools/S07604514_Padded.bed'

[bwa and picard]

[merge samples]

[remove sam reps]
[remove bam reps]

[mark duplicates]

[remove final sorted]

[base recalibrate]

[index dup]

[ApplyBQSR]

[remove dups]
