#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=139
#SBATCH --mem=2G
#SBATCH --job-name=plotseg%a
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/b1131/saya/bbcar/out/plot_cn_segments%a.out

cd /projects/b1131/saya/bbcar/

## Load GATK 
module purge all 
export PATH="/projects/b1131/saya/bbcar/tools/gatk-4.2.5.0:$PATH"
module load java/jdk1.8.0_191
module load R/3.3.1

## Directories etc.
DENOISED_CN_DIR="/projects/b1131/saya/bbcar/data/02b_cnv/03_denoised_counts"
ALLELIC_CTS_DIR="/projects/b1131/saya/bbcar/data/02b_cnv/05_allelic_counts"
CONTIGUOUS_CN_DIR="/projects/b1131/saya/bbcar/data/02b_cnv/06_contiguous_cn_segs"
CALLED_CNSEG_DIR="/projects/b1131/saya/bbcar/data/02b_cnv/07_called_cn_segs"
PLOT_DIR="/projects/b1131/saya/bbcar/plots/cnv/modeled_segments"

DICT=/projects/b1122/references/GRCH38/v0/Homo_sapiens_assembly38.dict

mkdir -p $PLOT_DIR
mkdir -p $PLOT_DIR/tissue_only/
mkdir -p $PLOT_DIR/tissue_normal/

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/bbcar/data/sample_ids_all_tissue.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

## Obtain IDs of samples processed at U Chicago (need to use different interval file) 
IFS=$'\n' read -d '' -r -a uchicago < /projects/b1131/saya/bbcar/data/sample_ids_uchicago.txt

## Obtain IDs of samples that have germline
IFS=$'\n' read -d '' -r -a germline < /projects/b1131/saya/bbcar/data/sample_ids_all_germline.txt

#### Run tissue-only for all samples ####
gatk PlotModeledSegments \
    --denoised-copy-ratios ${DENOISED_CN_DIR}/${sampleid}.denoisedCR.tsv \
    --allelic-counts ${CONTIGUOUS_CN_DIR}/tissue_only/${sampleid}.hets.tsv \
    --segments ${CONTIGUOUS_CN_DIR}/tissue_only/${sampleid}.modelFinal.seg \
    --sequence-dictionary ${DICT} \
    --minimum-contig-length 46709983 \
    --output ${PLOT_DIR}/tissue_only \
    --output-prefix ${sampleid}

#### Run for tissue-germline pairs if germline is present ####
if [[ " ${germline[*]} " =~ " ${sampleid} " ]]; then
    gatk PlotModeledSegments \
        --denoised-copy-ratios ${DENOISED_CN_DIR}/${sampleid}.denoisedCR.tsv \
        --allelic-counts ${CONTIGUOUS_CN_DIR}/tissue_normal/${sampleid}.hets.tsv \
        --segments ${CONTIGUOUS_CN_DIR}/tissue_normal/${sampleid}.modelFinal.seg \
        --sequence-dictionary ${DICT} \
        --minimum-contig-length 46709983 \
        --output ${PLOT_DIR}/tissue_normal \
        --output-prefix ${sampleid}
fi

