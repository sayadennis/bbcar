#!/bin/bash
#SBATCH -A p30791
#SBATCH -p short
#SBATCH -t 2:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=0-239
#SBATCH --mem=2G
#SBATCH --job-name=samqa
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/samtools_qa_%a.out

module purge all
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

## Set the input and output directories
aligndir='/projects/b1131/saya/new_bbcar/data/01_alignment'
metricsdir='/projects/b1131/saya/new_bbcar/metrics/samtools_qc'

mkdir -p $metricsdir

#### Run Samtools QC ####

if [[ " ${tissue[*]} " =~ " ${sampleid} " ]]; then
    echo "## Running tissue ##"
    samtools stats ${aligndir}/tissue/aligned/${sampleid}_bqsr.bam | grep ^SN | cut -f 2- > ${metricsdir}/${sampleid}_tissue_metrics.txt
fi

if [[ " ${germline[*]} " =~ " ${sampleid} " ]]; then
    echo "## Running germline ##"
    samtools stats ${aligndir}/germline/aligned/${sampleid}_bqsr.bam | grep ^SN | cut -f 2- > ${metricsdir}/${sampleid}_germline_metrics.txt
fi

