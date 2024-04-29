#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --array=0-239
#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH --mem=1G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="dpfilter"
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/DPfilter_%a.out

cd /projects/b1131/saya/new_bbcar/

## Load GATK 
module purge all
module load bcftools/1.10.1

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_tissue.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

## Obtain IDs of samples that have tissue or germline
IFS=$'\n' read -d '' -r -a tissue < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_tissue.txt
IFS=$'\n' read -d '' -r -a germline < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_germline.txt

function dp_filter() {
    local sampleid=$1
    local mode=$2

    if [ "$mode" = "germline_only" ]; then
        ext="haplotype"
    else
        ext="filtered"
    fi

    ## Set input and output directories
    din="/projects/b1131/saya/new_bbcar/data/02a_mutation/02_variant_calls/${mode}"
    dout=$din
    
    ## Apply filter on DP (depth)
    bcftools filter \
        --include "INFO/DP>=20" \
        --output-type v \
        --output $dout/${sampleid}_DPfiltered.vcf \
        $din/${sampleid}_${ext}.vcf
    
    ## If tissue sample, apply filter to correct for FFPE artifacts
    ## Classical filter taken from PMID 34261476 (Winham et al)
    if [ "$mode" != "germline_only" ]; then
        bcftools filter \
            --exclude '(FORMAT/AF<0.05) | (FORMAT/AF<0.10 && REF="C" && ALT="T")' \
            --output-type v \
            --output $dout/${sampleid}_DPfiltered_classicalAF.vcf \
            $din/${sampleid}_DP${ext}.vcf
    fi
}

## Tissue only
if [[ " ${tissue[*]} " =~ " ${sampleid} " ]]; then
    echo "########## Filtering variants on tissue sample ${sampleid} ##########"
    dp_filter ${sampleid} 'tissue_only'
else
    echo "########## Patient ID ${sampleid} is not in the tissue sample list ##########"
fi

## Tissue-normal and germline-only
if [[ " ${germline[*]} " =~ " ${sampleid} " ]]; then
    echo "########## Filtering variants on tissue-germline pair of sample ${sampleid} ##########"
    dp_filter ${sampleid} 'tissue_normal'
    echo "########## Filtering variants on germline sample ${sampleid} ##########"
    dp_filter ${sampleid} 'germline_only'
else
    echo "########## Patient ID ${sampleid} is not in the germline sample list ##########"
fi

