#!/bin/bash
#SBATCH -A p30791
#SBATCH -p normal
#SBATCH --array=0-239
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --mem=3G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="mutect2_%a"
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/call_variants_%a.out

cd /projects/b1131/saya/new_bbcar/

## Set GATK command with singularity image
module purge all
module load singularity

gatk() {
    singularity exec -B /projects:/projects /projects/p30791/gatk_4.5.0.sif gatk "$@"
}

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_tissue.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

## Obtain IDs of samples that have tissue or germline
IFS=$'\n' read -d '' -r -a tissue < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_tissue.txt
IFS=$'\n' read -d '' -r -a germline < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_germline.txt

## Obtain IDs of samples processed at U Chicago (need to use different interval file) 
IFS=$'\n' read -d '' -r -a uchicago < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_batch_1.txt

## Set reference, interval, germline resource directories
ref='/projects/p30791/hg38_ref/hg38.fa'
germres='/projects/b1131/saya/bbcar/genome_resources/GATK/af-only-gnomad.hg38.vcf.gz'
int_dir='/projects/b1131/saya/bbcar/interval_lists'

## Instructions taken from https://gatk.broadinstitute.org/hc/en-us/articles/360035531132 
## Under section "A step-by-step guide to the new Mutect2 Read Orientation Artifacts Workflow"

# set PON file name
pon='/projects/b1131/saya/bbcar/data/02a_mutation/02_variant_calls/germline_only/bbcar_pon.vcf.gz'

function call_variants() {
    local sampleid=$1  # 1004
    local mode=$2  # tissue_only, tissue_normal, or germline_only
    local interval=$3  # interval filepath

    if [ "$mode" = "germline_only" ]; then
        din='/projects/b1131/saya/new_bbcar/data/01_alignment/germline/aligned'
        dout='/projects/b1131/saya/new_bbcar/data/02a_mutation/02_variant_calls/germline_only'

        mkdir -p $dout

        ## Create output with raw data used to learn the orientation bias model
        gatk --java-options "-Xmx3g" HaplotypeCaller \
                -R $ref \
                --intervals $interval \
                -I $din/${sampleid}_bqsr.bam \
                -O $dout/${sampleid}_haplotype.vcf

    else
        if [ "$mode" = "tissue_normal" ]; then
            ## Set input and output directories 
            din_tiss='/projects/b1131/saya/new_bbcar/data/01_alignment/tissue/aligned'
            din_germ='/projects/b1131/saya/new_bbcar/data/01_alignment/germline/aligned'
            dout='/projects/b1131/saya/new_bbcar/data/02a_mutation/02_variant_calls/tissue_normal'

            mkdir -p $dout/interim/

            ## Create output with raw data used to learn the orientation bias model
            gatk Mutect2 \
                 -R $ref \
                 -L $interval \
                 -I $din_tiss/${sampleid}_bqsr.bam \
                 -I $din_germ/${sampleid}_bqsr.bam \
                 --normal-sample ${sampleid}_germline \
                 -germline-resource $germres \
                 -pon $pon \
                 --f1r2-tar-gz $dout/interim/${sampleid}_f1r2.tar.gz \
                 -O $dout/interim/${sampleid}_unfiltered.vcf
        elif [ "$mode" = "tissue_only" ]; then
            ## Set input and output directories 
            din_tiss='/projects/b1131/saya/new_bbcar/data/01_alignment/tissue/aligned'
            din_germ=''
            dout='/projects/b1131/saya/new_bbcar/data/02a_mutation/02_variant_calls/tissue_only'

            mkdir -p $dout/interim

            ## Create output with raw data used to learn the orientation bias model
            gatk Mutect2 \
                -R $ref \
                -L $interval \
                -I $din_tiss/${sampleid}_bqsr.bam \
                -germline-resource $germres \
                -pon $pon \
                --f1r2-tar-gz $dout/interim/${sampleid}_f1r2.tar.gz \
                -O $dout/interim/${sampleid}_unfiltered.vcf
        fi

        ## Pass this raw data to LearnReadOrientationModel
        gatk LearnReadOrientationModel \
            -I $dout/interim/${sampleid}_f1r2.tar.gz \
            -O $dout/interim/${sampleid}_read-orientation-model.tar.gz
        
        ## Run GetPileupSummaries to summarize read support for a set number of known variant sites.
        gatk GetPileupSummaries \
            --input $din_tiss/${sampleid}_bqsr.bam \
            --variant $germres \
            --intervals $interval \
            --output $dout/interim/${sampleid}_getpileupsummaries.table
        
        ## Estimate contamination with CalculateContamination.
        gatk CalculateContamination \
            -I $dout/interim/${sampleid}_getpileupsummaries.table \
            -tumor-segmentation $dout/interim/${sampleid}_segments.table \
            -O $dout/interim/${sampleid}_calculatecontamination.table
        
        ## Finally, pass the learned read orientation model to FilterMutectCallswith the -ob-priors argument:
        gatk FilterMutectCalls \
            -V $dout/interim/${sampleid}_unfiltered.vcf \
            -R $ref \
            --ob-priors $dout/interim/${sampleid}_read-orientation-model.tar.gz \
            -O $dout/${sampleid}_filtered.vcf
    fi
}

#### Run ####

if [[ " ${uchicago[*]} " =~ " ${sampleid} " ]]; then
    intfile=${int_dir}/SureSelect_v5/hg38.v5.preprocessed.interval_list
else
    intfile=${int_dir}/SureSelect_v6/hg38.v6.preprocessed.interval_list
fi

## Tissue only
if [[ " ${tissue[*]} " =~ " ${sampleid} " ]]; then
    echo "########## Calling variants on tissue sample ${sampleid} ##########"
    call_variants ${sampleid} 'tissue_only' $intfile
else
    echo '########## Patient ID ${sampleid} is not in the tissue sample list ##########'
fi

## Tissue-normal and germline-only
if [[ " ${germline[*]} " =~ " ${sampleid} " ]]; then
    echo "########## Calling variants on tissue-germline pair of sample ${sampleid} ##########"
    call_variants ${sampleid} 'tissue_normal' $intfile
    echo "########## Calling variants on germline sample ${sampleid} ##########"
    call_variants ${sampleid} 'germline_only' $intfile
else
    echo '########## Patient ID ${sampleid} is not in the germline sample list ##########'
fi

