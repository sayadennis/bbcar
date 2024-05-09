#!/bin/bash
#SBATCH -A p31931
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 5
#SBATCH --array=0-239
#SBATCH --mem=50G
#SBATCH --job-name=bwa
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/align_bwa_%a.out

module purge all
module load bwa
module load picard/2.6.0
module load gatk/4.1.0
module load samtools

cd /projects/b1131/saya/new_bbcar/

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_all.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

## Obtain IDs of samples that have tissue or germline
IFS=$'\n' read -d '' -r -a tissue < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_tissue.txt
IFS=$'\n' read -d '' -r -a germline < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_germline.txt

## Obtain IDs of samples processed at U Chicago (need to use different interval file) 
IFS=$'\n' read -d '' -r -a batch1 < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_batch_1.txt

set -u -e -x

## Reference genome
FA='/projects/p30791/hg38_ref/hg38.fa'

# index the reference genome
# bwa index -p /projects/p30791/hg38_ref/hg38.fa $FA

## Software and interval locations
pic='java -jar /software/picard/2.6.0/picard-tools-2.6.0/picard.jar'
dbsnp='/projects/b1131/saya/bbcar/genome_resources/Homo_sapiens_assembly38.dbsnp138.vcf' # is this file appropriate for this?
gold1000Indel='/projects/b1131/saya/bbcar/genome_resources/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz' # is this file appropriate for this?

## If tissue comes from batch 1, then set tissue intervals for V5; otherwise set to V6
if [[ " ${batch1[*]} " =~ " ${sampleid} " ]]; then
    tissue_interval='/projects/b1122/gannon/bbcar/RAW_data/int_lst/SureSelect_v5/hg38/hg38.preprocessed.interval_list' # V5
    germline_interval='/projects/b1122/gannon/bbcar/RAW_data/int_lst/SureSelect_v6/hg38.preprocessed.interval_list' # V6
else
    tissue_interval='/projects/b1122/gannon/bbcar/RAW_data/int_lst/SureSelect_v6/hg38.preprocessed.interval_list' # V6
    germline_interval='/projects/b1122/gannon/bbcar/RAW_data/int_lst/SureSelect_v6/hg38.preprocessed.interval_list' # V6
fi

## Set the input and output directories
fastqdir='/projects/b1131/saya/new_bbcar/data/00_trimmed'
aligndir='/projects/b1131/saya/new_bbcar/data/01_alignment'

mkdir -p $aligndir
for tissuetype in tissue germline; do
    for dirname in metrics rec_tables interim aligned; do
        mkdir -p "${aligndir}/${tissuetype}/${dirname}"
    done
done

function run_alignment() {
    local sampleid=$1  # 1004
    local tissuetype=$2  # tissue -or- germline
    local interval=$3  # interval filepath

    # Define output directories
    metrics=${aligndir}/${tissuetype}/metrics
    rec_tables=${aligndir}/${tissuetype}/rec_tables
    interim=${aligndir}/${tissuetype}/interim
    aligned=${aligndir}/${tissuetype}/aligned

    # First, align the paired FASTQ reads
    for fpattern in $(ls ${fastqdir}/${tissuetype}/${sampleid}_*_paired.fastq.gz \
        | xargs -n 1 basename \
        | awk -F'[_.]' '{print $1"_"$2"_"$3}' \
        | sed 's/_R[12]$//' \
        | sort -u);
    do  
        readarray -t files_array <<< "$(ls ${fastqdir}/${tissuetype}/${fpattern}_*_paired.fastq.gz)"
    
        if [[ ${#files_array[@]} -gt 2 ]]; then
            echo "Error: there are more than 2 files that match file pattern: ${fastqdir}/${tissuetype}/${fpattern}_*_paired.fastq.gz"
            exit 1
        fi  
    
        R1=${files_array[0]}
        R2=${files_array[1]}

        header=$(zcat $R1 | head -n 1)
        id=$(echo $header | head -n 1 | cut -f 3-4 -d':' | sed 's/@//' | sed 's/:/_/g')
        echo "@RG\tID:${id}\tSM:${sampleid}_${tissuetype}\tLB:library1\tPL:ILLUMINA\tPU:${id}"
    
        bwa mem -M -t 24 \
            -R $(echo "@RG\tID:${id}\tSM:${sampleid}_${tissuetype}\tLB:library1\tPL:ILLUMINA\tPU:${id}") \
            $FA $R1 $R2 > $interim/${fpattern}.sam
    
        $pic SortSam I=${interim}/${fpattern}.sam O=${interim}/${fpattern}_sorted.bam SORT_ORDER=coordinate
    done

    # Merge the BAM files if there are multiple
    $pic MergeSamFiles CREATE_INDEX=true USE_THREADING=true \
        $(for x in $(ls -1 ${interim}/${sampleid}_*_sorted.bam); do echo -n "I=$x "; done) \
        O=${interim}/${sampleid}_final_sorted.bam
    
    for fpattern in $(ls ${fastqdir}/${tissuetype}/${sampleid}_*_paired.fastq.gz \
        | xargs -n 1 basename \
        | awk -F'[_.]' '{print $1"_"$2"_"$3}' \
        | sed 's/_R[12]$//' \
        | sort -u);
    do  
        rm -f ${interim}/${fpattern}.sam
        rm -f ${interim}/${fpattern}_sorted.bam
    done
    
    # Mark duplicated in the merged SAM files
    $pic MarkDuplicates \
        I=${interim}/${sampleid}_final_sorted.bam \
        O=${interim}/${sampleid}_dup.bam \
        M=${metrics}/${sampleid}_${tissuetype}_reads.mdup.metrics.txt
    
    rm -f ${interim}/${sampleid}_final_sorted.bam
    rm -f ${interim}/${sampleid}_final_sorted.bai
    
    # Base recalibration
    gatk BaseRecalibrator \
        -I ${interim}/${sampleid}_dup.bam \
        -R $FA \
        --known-sites $dbsnp \
        --known-sites $gold1000Indel \
        -O ${rec_tables}/${sampleid}_${tissuetype}_recal_data.table
    
    # Index bam files
    samtools index ${interim}/${sampleid}_dup.bam
    
    # Apply BQSR
    gatk ApplyBQSR \
        -R $FA \
        -I ${interim}/${sampleid}_dup.bam \
        --bqsr-recal-file ${rec_tables}/${sampleid}_${tissuetype}_recal_data.table \
        -L ${interval} \
        -O ${aligned}/${sampleid}_bqsr.bam
    
    rm -f ${interim}/${sampleid}_dup.bam
    rm -f ${interim}/${sampleid}_dup.bam.bai
}

#######################
#### Run alignment ####
#######################

if [[ " ${tissue[*]} " =~ " ${sampleid} " ]]; then
    echo "########## Running alignment on tissue sample reads for patient ID ${sampleid} ##########"
    run_alignment ${sampleid} 'tissue' $tissue_interval
else
    echo "########## Patient ID ${sampleid} is not in the tissue sample list ##########"
fi

if [[ " ${germline[*]} " =~ " ${sampleid} " ]]; then
    echo "########## Running alignment on germline sample reads for patient ID ${sampleid} ##########"
    run_alignment ${sampleid} 'germline' $germline_interval
else
    echo "########## Patient ID ${sampleid} is not in the germline sample list ##########"
fi

