#!/bin/bash
#SBATCH -A p31931
#SBATCH -p short
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH --array=0-239
#SBATCH --mem=2G
#SBATCH --job-name=alleleCount
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/run_alleleCounter%a.out

module purge all
module use ${HOME}/mymodules/
module load alleleCounter

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_all.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

## Obtain IDs of samples that have germline
IFS=$'\n' read -d '' -r -a tissue < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_tissue.txt
IFS=$'\n' read -d '' -r -a germline < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_germline.txt

## Obtain IDs of samples processed at U Chicago (need to use different interval file) 
IFS=$'\n' read -d '' -r -a uchicago < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_batch_1.txt

outdir=/projects/b1131/saya/new_bbcar/data/02b_cnv/signatures/01_alleleCounts
mkdir -p $outdir

ref_fasta=/projects/p30791/hg38_ref/hg38.fa

tissue_seqfile=/projects/b1131/saya/new_bbcar/data/01_alignment/tissue/aligned/${sampleid}_bqsr.bam
germline_seqfile=/projects/b1131/saya/new_bbcar/data/01_alignment/germline/aligned/${sampleid}_bqsr.bam

## Select the correct interval based on samples sequenced at U of Chicago (TREAT GERMLINE DIFFERENTLY?! - CHECK ALIGNMENT SCRIPT)
if [[ " ${uchicago[*]} " =~ " ${sampleid} " ]]; then
    int_file='/projects/b1122/gannon/bbcar/RAW_data/int_lst/SureSelect_v5/hg38/hg38.preprocessed.interval_list' # V5
    loci_file=${outdir}/v5interval_loci_format_for_alleleCounter.bed
else
    int_file='/projects/b1122/gannon/bbcar/RAW_data/int_lst/SureSelect_v6/hg38.preprocessed.interval_list' # V6
    loci_file=${outdir}/v6interval_loci_format_for_alleleCounter.bed
fi

## Create loci file for alleleCounter if it doesn't exist already 
if [ ! -f $loci_file ]
then
    echo "Loci file does not exist. Creating $loci_file from $int_file"
    grep "^chr" $int_file > $loci_file
else
    echo "Loci file already exists. Using $loci_file"
fi

## Run alleleCounter for tissue and germline 
if [[ " ${germline[*]} " =~ " ${sampleid} " ]] && [[ " ${tissue[*]} " =~ " ${sampleid} " ]]; then
    echo "Running alleleCounter for tissue sequences $tissue_seqfile"
    alleleCounter \
        -b $tissue_seqfile \
        -l $loci_file \
        -o ${outdir}/${sampleid}_tissue_alleleFrequencies.txt \
        -m 20 \
        -q 35 \
        --dense-snps \
        -r $ref_fasta

    echo "Running alleleCounter for germline sequences $germline_seqfile"
    alleleCounter \
        -b $germline_seqfile \
        -l $loci_file \
        -o ${outdir}/${sampleid}_germline_alleleFrequencies.txt \
        -m 20 \
        -q 35 \
        --dense-snps \
        -r $ref_fasta
fi
