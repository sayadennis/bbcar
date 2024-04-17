#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=0-201
#SBATCH --mem=12G
#SBATCH --job-name=callseg%a
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/b1131/saya/bbcar/out/call_cn_segments%a.out

cd /projects/b1131/saya/bbcar/

## Load GATK 
module purge all 
export PATH="/projects/b1131/saya/bbcar/tools/gatk-4.2.5.0:$PATH"
module load java/jdk1.8.0_191

## Directories etc.
CONTIGUOUS_CN_DIR="/projects/b1131/saya/bbcar/data/02b_cnv/06_contiguous_cn_segs"
CALLED_CNSEG_DIR="/projects/b1131/saya/bbcar/data/02b_cnv/07_called_cn_segs"

mkdir -p $CALLED_CNSEG_DIR
mkdir -p $CALLED_CNSEG_DIR/tissue_only/
mkdir -p $CALLED_CNSEG_DIR/tissue_normal/

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/bbcar/data/sample_ids_all_tissue.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

## Obtain IDs of samples processed at U Chicago (need to use different interval file) 
IFS=$'\n' read -d '' -r -a uchicago < /projects/b1131/saya/bbcar/data/sample_ids_uchicago.txt

## Obtain IDs of samples that have germline
IFS=$'\n' read -d '' -r -a germline < /projects/b1131/saya/bbcar/data/sample_ids_all_germline.txt

#### Run tissue-only for all samples ####
gatk CallCopyRatioSegments \
    --input ${CONTIGUOUS_CN_DIR}/tissue_only/${sampleid}.cr.seg \
    --output ${CALLED_CNSEG_DIR}/tissue_only/${sampleid}.called.seg

#### Run for tissue-germline pairs if germline is present ####
if [[ " ${germline[*]} " =~ " ${sampleid} " ]]; then
    gatk CallCopyRatioSegments \
        --input ${CONTIGUOUS_CN_DIR}/tissue_normal/${sampleid}.cr.seg \
        --output ${CALLED_CNSEG_DIR}/tissue_normal/${sampleid}.called.seg
fi

