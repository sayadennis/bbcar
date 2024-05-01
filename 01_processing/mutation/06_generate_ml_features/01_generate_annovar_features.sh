#!/bin/bash
#SBATCH -A p30791
#SBATCH -p short
#SBATCH --array=0-239
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 3:00:00
#SBATCH --mem=1G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=avfts_%a
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/generate_annovar_features_%a.out

module purge all
module load python-miniconda3/4.12.0
source activate bbcarenv

cd ${HOME}/bbcar/repo/01_processing/mutation/06_generate_ml_features/

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_tissue.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

## Obtain IDs of samples that have germline
IFS=$'\n' read -d '' -r -a tissue < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_tissue.txt
IFS=$'\n' read -d '' -r -a germline < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_germline.txt

## Tissue only
if [[ " ${tissue[*]} " =~ " ${sampleid} " ]]; then
    python 01_generate_annovar_features.py $sampleid "tissue_only"
else
    echo '########## Patient ID ${sampleid} is not in the tissue sample list ##########'
fi

## Tissue-normal and germline-only
if [[ " ${germline[*]} " =~ " ${sampleid} " ]]; then
    python 01_generate_annovar_features.py $sampleid "tissue_normal"
    python 01_generate_annovar_features.py $sampleid "germline_only"
else
    echo '########## Patient ID ${sampleid} is not in the germline sample list ##########'
fi

