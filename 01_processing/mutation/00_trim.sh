#!/bin/bash
#SBATCH -A p31931
#SBATCH -p short
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --array=0-239
#SBATCH --mem=16G
#SBATCH --job-name=trim
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/trim_%a.out

cd /projects/b1131/saya/new_bbcar/

module purge all
module load java/jdk1.8.0_191
module load trimmomatic/0.39

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_all.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

## Obtain IDs of samples that have tissue or germline
IFS=$'\n' read -d '' -r -a tissue < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_tissue.txt
IFS=$'\n' read -d '' -r -a germline < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_germline.txt

## Obtain IDs of samples processed at U Chicago (need to use different interval file) 
IFS=$'\n' read -d '' -r -a batch3 < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_batch_3.txt

## Set the input and output directories
fastqdir='/projects/b1131/saya/new_bbcar/data/fastq'
trimdir='/projects/b1131/saya/new_bbcar/data/00_trimmed'
adapterfile='/projects/b1131/saya/bbcar/genome_resources/TruSeq3-PE-2.fa'

mkdir -p ${trimdir}/tissue/
mkdir -p ${trimdir}/germline/

function trim_reads() {
    local sampleid=$1  # 1004
    local tissuetype=$2  # tissue -or- germline

    for fpattern in $(ls ${fastqdir}/${tissuetype}/${sampleid}_*.fastq.gz \
        | xargs -n 1 basename \
        | awk -F'[_.]' '{print $1"_"$2"_"$3}' \
        | sed 's/_R[12]$//' \
        | sort -u);
    do  
        readarray -t files_array <<< "$(ls ${fastqdir}/${tissuetype}/${fpattern}*)"
    
        if [[ ${#files_array[@]} -gt 2 ]]; then
            echo "Error: there are more than 2 files that match file pattern: ${fastqdir}/${tissuetype}/${fpattern}"
            exit 1
        else
            echo "Running Trimmomatic for files with name pattern: ${fastqdir}/${tissuetype}/${fpattern}*"
        fi  
    
        if [[ "${batch3[*]} " =~ " ${sampleid} " ]]; then
            cropcommand="HEADCROP:6"
            echo "Patient ID is in batch 3 - cropping the leading 6bp."
        else
            cropcommand="CROP:72"
            echo "Patient ID is not in batch 3 - cropping the trailing 2bp."
        fi

        R1=${files_array[0]}
        R2=${files_array[1]}
    
        java -jar /software/trimmomatic/0.39/trimmomatic-0.39.jar \
            PE \
            -threads 8 \
            $R1 $R2 \
            ${trimdir}/${tissuetype}/${fpattern}_R1_trimmed_paired.fastq.gz ${trimdir}/${tissuetype}/${fpattern}_R1_trimmed_unpaired.fastq.gz \
            ${trimdir}/${tissuetype}/${fpattern}_R2_trimmed_paired.fastq.gz ${trimdir}/${tissuetype}/${fpattern}_R2_trimmed_unpaired.fastq.gz \
            ILLUMINACLIP:${adapterfile}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 $(echo ${cropcommand})
    done
}

if [[ " ${tissue[*]} " =~ " ${sampleid} " ]]; then
    echo "########## Trimming reads for tissue sample for patient ID ${sampleid} ##########"
    trim_reads ${sampleid} 'tissue'
else
    echo "########## Patient ID ${sampleid} is not in the tissue sample list ##########"
fi

if [[ " ${germline[*]} " =~ " ${sampleid} " ]]; then
    echo "########## Trimming reads for germline sample for patient ID ${sampleid} ##########"
    trim_reads ${sampleid} 'germline'
else
    echo "########## Patient ID ${sampleid} is not in the germline sample list ##########"
fi

