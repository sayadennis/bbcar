#!/bin/bash                                                                                                  
[settings]

module purge all
module load samtools
module load java/jdk1.8.0_25
module load gatk/4.1.0
module load R
module load blas-lapack/3.5.0

cd /projects/b1042/ClareLab/GATK_PON/PON/germ_readcount

refDir=/projects/p30007/Zexian/reference/hg19/ucsc.hg19.fasta
dictDir=/projects/p30007/Zexian/reference/hg19/ucsc.hg19.dict
idxDir=/projects/p30007/Zexian/reference/hg19/ucsc.hg19.fasta.fai
int_lst=/projects/p30007/gannon/bbcar/new_int
germ_bam=/projects/b1042/ClareLab/GATK_PON/bam/germ
tissue_bam=/projects/b1042/ClareLab/GATK_PON/bam/tissue
processed_int_lst=/projects/p30007/gannon/bbcar/gatk_out
PONDir=/projects/b1042/ClareLab/GATK_PON/PON
segment=/projects/b1042/ClareLab/GATK_PON/PON/segs
germOut=/projects/b1042/ClareLab/GATK_PON/PON/germ_readcount
tissOut=/projects/b1042/ClareLab/GATK_PON/PON/tiss_readcount
denoiseOut=/projects/b1042/ClareLab/GATK_PON/denoised_cnvs
plots=/projects/b1042/ClareLab/GATK_PON/plots

                                              
gatk --java-options "-Xmx12g" DenoiseReadCounts -I $tissOut/[patid].counts.hdf5 --count-panel-of-normals $PONDir/germtest.pon.hdf5 --standardized-copy-ratios $denoiseOut/[patid]_clean.standardizedCR.tsv --denoised-copy-ratios $denoiseOut/[patid]_clean.denoisedCR.tsv

mkdir $plots/[patid]
gatk PlotDenoisedCopyRatios --standardized-copy-ratios $denoiseOut/[patid]_clean.standardizedCR.tsv --denoised-copy-ratios $denoiseOut/[patid]_clean.denoisedCR.tsv --sequence-dictionary $dictDir --minimum-contig-length 48129895 --output $plots/[patid] --output-prefix [patid]_clean


gatk --java-options "-Xmx100g" CollectAllelicCounts -L $processed_int_lst/S07604514_padded.targets.interval_list -I $tissue_bam/[patid]_bqsr.bam -R $refDir -O $PONDir/refCounts/[patid]_clean.allelicCounts.tsv

mkdir $segment/[patid]
gatk --java-options "-Xmx50g" ModelSegments --denoised-copy-ratios $denoiseOut/[patid]_clean.denoisedCR.tsv --output $segment/[patid]  --output-prefix [patid]_clean

gatk CallCopyRatioSegments --input $segment/[patid]/[patid]_clean.cr.seg  --output $segment/called/[patid]_clean.called.seg


gatk PlotModeledSegments --denoised-copy-ratios $denoiseOut/[patid]_clean.denoisedCR.tsv --segments $segment/[patid]/[patid]_clean.modelFinal.seg --sequence-dictionary $dictDir --minimum-contig-length 48129895  --output $plots/[patid]  --output-prefix [patid]_clean
