#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --array=0-201
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 4:00:00
#SBATCH --mem=10G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=av_to%a
#SBATCH --output=/projects/b1131/saya/bbcar/out/annovar_tumor_only_%a.out

cd /projects/b1131/saya/bbcar/

## Load necessary modules 
module purge all
module load perl/5.16

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/bbcar/data/02a_mutation/sample_names_tumor_only.txt

## Set input and output directories 
din='/projects/b1131/saya/bbcar/data/02a_mutation/02_variant_calls/tumor_only'
dav='/projects/b1131/saya/bbcar/tools/annovar'

for filter_type in liberal classical strict; do
    dout=/projects/b1131/saya/bbcar/exploratory_winham_filter/${filter_type}/03_annotated_variants/annovar/tumor_only
    mkdir -p $dout
    pon_source='bbcar'
    fin=${input_args[$SLURM_ARRAY_TASK_ID]}_DPfiltered_${pon_source}pon.vcf
    fout=${input_args[$SLURM_ARRAY_TASK_ID]}_${pon_source}pon

    perl ${dav}/table_annovar.pl \
        ${din}/${fin} \
        ${dav}/humandb/ \
        -buildver hg38 \
        -out $dout/$fout \
        -remove \
        -protocol refGene,knownGene,ensGene,avsnp150,dbnsfp35a,dbnsfp31a_interpro,exac03,gnomad211_exome \
        -operation g,g,g,f,f,f,f,f \
        -nastring . -vcfinput
    #
done
