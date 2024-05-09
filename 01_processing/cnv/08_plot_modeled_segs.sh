#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -t 1:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --array=0-239
#SBATCH --mem=2G
#SBATCH --job-name=plotseg%a
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/plot_cn_segments%a.out

cd /projects/b1131/saya/new_bbcar/

## Load GATK 
module purge all 
module load R/3.3.1
module load singularity

gatk() {
    singularity exec -B /projects:/projects /projects/p30791/gatk_4.5.0.sif gatk "$@"
}

## Directories etc.
DENOISED_CN_DIR="/projects/b1131/saya/new_bbcar/data/02b_cnv/03_denoised_counts"
ALLELIC_CTS_DIR="/projects/b1131/saya/new_bbcar/data/02b_cnv/05_allelic_counts"
CONTIGUOUS_CN_DIR="/projects/b1131/saya/new_bbcar/data/02b_cnv/06_contiguous_cn_segs"
CALLED_CNSEG_DIR="/projects/b1131/saya/new_bbcar/data/02b_cnv/07_called_cn_segs"
PLOT_DIR="/projects/b1131/saya/new_bbcar/plots/cnv/modeled_segments"

DICT=/projects/b1122/references/GRCH38/v0/Homo_sapiens_assembly38.dict

mkdir -p $PLOT_DIR/tissue_only/
mkdir -p $PLOT_DIR/tissue_normal/

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_all.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

## Obtain IDs of samples that have germline
IFS=$'\n' read -d '' -r -a tissue < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_tissue.txt
IFS=$'\n' read -d '' -r -a germline < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_germline.txt

## Obtain IDs of samples processed at U Chicago (need to use different interval file) 
IFS=$'\n' read -d '' -r -a uchicago < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_batch_1.txt

#### Run tissue-only for all samples ####
if [[ " ${tissue[*]} " =~ " ${sampleid} " ]]; then
    gatk PlotModeledSegments \
        --denoised-copy-ratios ${DENOISED_CN_DIR}/${sampleid}.denoisedCR.tsv \
        --allelic-counts ${CONTIGUOUS_CN_DIR}/tissue_only/${sampleid}.hets.tsv \
        --segments ${CONTIGUOUS_CN_DIR}/tissue_only/${sampleid}.modelFinal.seg \
        --sequence-dictionary ${DICT} \
        --minimum-contig-length 46709983 \
        --output ${PLOT_DIR}/tissue_only \
        --output-prefix ${sampleid}
fi

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

