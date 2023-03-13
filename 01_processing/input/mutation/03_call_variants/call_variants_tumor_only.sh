#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --array=0-201
#SBATCH -n 1
#SBATCH -t 4:00:00
#SBATCH --mem=1G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name="mutect2_to_%a"
#SBATCH --output=/projects/b1131/saya/bbcar/out/call_variants_tumor_only_%a.out

cd /projects/b1131/saya/bbcar/

## Load GATK 
module purge all
export PATH="/projects/b1131/saya/bbcar/tools/gatk-4.2.5.0:$PATH"

#### Create sample names TXT file for job array (RUN THIS IN ADVANCE) ####
# ls /projects/b1131/saya/bbcar/01_alignment/tissue/aligned/*.bam | xargs -n1 basename | tr "_" "\n" | grep "t" > /projects/b1131/saya/bbcar/sample_names_tumor_only.txt

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/bbcar/data/02a_mutation/sample_names_tumor_only.txt

## Set input and output directories
din='/projects/b1131/saya/bbcar/data/01_alignment/tissue/aligned'
dout='/projects/b1131/saya/bbcar/data/02a_mutation/02_variant_calls/tumor_only'

## Set reference, interval, germline resource, and PON filenames 
ref='/projects/p30791/hg38_ref/hg38.fa'
interval='/projects/b1122/gannon/bbcar/RAW_data/int_lst/SureSelect_v6/hg38.preprocessed.interval_list'
germres='/projects/b1131/saya/bbcar/genome_resources/GATK/af-only-gnomad.hg38.vcf.gz'

## Instructions taken from https://gatk.broadinstitute.org/hc/en-us/articles/360035531132 
## Under section "A step-by-step guide to the new Mutect2 Read Orientation Artifacts Workflow"

for pon_source in bbcar 1000g; do
        # set PON file name
        if [ "$pon_source" = "bbcar" ]
        then pon='/projects/b1131/saya/bbcar/02_variant_calls/germline_only/bbcar_pon.vcf.gz'
        elif [ "$pon_source" = "1000g" ]
        then pon='/projects/b1131/saya/bbcar/genome_resources/GATK/1000g_pon.hg38.vcf.gz'
        fi
        # 

        ## Create output with raw data used to learn the orientation bias model
        gatk Mutect2 \
                -R $ref \
                -L $interval \
                -I $din/${input_args[$SLURM_ARRAY_TASK_ID]}_bqsr.bam \
                -germline-resource $germres \
                -pon $pon \
                --f1r2-tar-gz $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_f1r2_${pon_source}pon.tar.gz \
                -O $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_unfiltered_${pon_source}pon.vcf
        #

        ## Pass this raw data to LearnReadOrientationModel
        gatk LearnReadOrientationModel -I $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_f1r2_${pon_source}pon.tar.gz -O $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_read-orientation-model_${pon_source}pon.tar.gz

        ## Run GetPileupSummaries to summarize read support for a set number of known variant sites.
        gatk GetPileupSummaries \
        --input $din/${input_args[$SLURM_ARRAY_TASK_ID]}_bqsr.bam \
        --variant $germres \
        --intervals $interval \
        --output $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_getpileupsummaries_${pon_source}pon.table
        #

        ## Estimate contamination with CalculateContamination.
        gatk CalculateContamination \
                -I $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_getpileupsummaries_${pon_source}pon.table \
                -tumor-segmentation $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_segments_${pon_source}pon.table \
                -O $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_calculatecontamination_${pon_source}pon.table
        #

        ## Finally, pass the learned read orientation model to FilterMutectCallswith the -ob-priors argument:
        gatk FilterMutectCalls -V $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_unfiltered_${pon_source}pon.vcf \
                -R $ref \
                --ob-priors $dout/interim/${input_args[$SLURM_ARRAY_TASK_ID]}_read-orientation-model_${pon_source}pon.tar.gz \
                -O $dout/${input_args[$SLURM_ARRAY_TASK_ID]}_filtered_${pon_source}pon.vcf
                # [--tumor-segmentation $dout/interim/segments.table] \
                # [--contamination-table contamination.table] \
        #
done
