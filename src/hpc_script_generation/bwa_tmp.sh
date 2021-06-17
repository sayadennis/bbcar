#!/bin/bash              
[settings]

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
