#!/bin/bash
#SBATCH -A p30791
#SBATCH -p short
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --array=0-239
#SBATCH --mem=50G
#SBATCH --job-name=fastqc
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/fastqc_%a.out

module purge all
module load java/jdk1.8.0_191
module load openjdk
module load fastqc/0.12.0
module load multiqc/1.2

indir="/projects/b1131/saya/new_bbcar/data/00_trimmed"
outdir="/projects/b1131/saya/new_bbcar/metrics/fastqc_trimmed"
scratchdir="/scratch/srd6051/fastqc_trimmed"

cd $indir

mkdir -p $outdir
mkdir -p $scratchdir

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_all.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

## Obtain IDs of samples that have tissue or germline
IFS=$'\n' read -d '' -r -a tissue < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_tissue.txt
IFS=$'\n' read -d '' -r -a germline < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_germline.txt

if [[ " ${tissue[*]} " =~ " ${sampleid} " ]]; then
    fastqc \
        --outdir $outdir \
        --noextract \
        --memory 10000 \
        --threads 10 \
        --java /software/java/jdk1.8.0_191/bin/java \
        --dir $scratchdir \
        $(for x in $(ls -1 ./tissue/${sampleid}_*_trimmed_paired.fastq.gz); do echo -n "$x "; done)
fi

if [[ " ${germline[*]} " =~ " ${sampleid} " ]]; then
    fastqc \
        --outdir $outdir \
        --noextract \
        --memory 10000 \
        --threads 10 \
        --java /software/java/jdk1.8.0_191/bin/java \
        --dir $scratchdir \
        $(for x in $(ls -1 ./germline/${sampleid}_*_trimmed_paired.fastq.gz); do echo -n "$x "; done)
fi

