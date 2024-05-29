#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 48:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=1G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="somvcf"
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/create_somatic_vcfs.out

module purge all
module load python-miniconda3/4.12.0
source activate bbcarenv

cd /projects/b1131/saya/new_bbcar/

## Obtain IDs of samples that have tissue or germline
IFS=$'\n' read -d '' -r -a tissue < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_tissue.txt
IFS=$'\n' read -d '' -r -a germline < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_germline.txt

#### Copy the tissue-normal matched calls (100% somatic) ####
din="/projects/b1131/saya/new_bbcar/data/02a_mutation/02_variant_calls"
dout="/projects/b1131/saya/new_bbcar/data/02a_mutation/07_predicted_somatic/vcfs"

mkdir -p ${dout}

cp ${din}/tissue_normal/*_DPfiltered_classicalAF.vcf ${dout}/

#### Filter the tissue-only calls and select predicted-somatic ####
python ~/bbcar/repo/01_processing/mutation/07_feature_matrix/01_create_somatic_vcfs.py

