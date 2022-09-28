#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --array=0-51
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 4:00:00
#SBATCH --mem=20G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=av_tn_${SLURM_ARRAY_TASK_ID}
#SBATCH --output=bbcar/out/annovar_tumor_normal_%a.out

cd /projects/b1131/saya/bbcar/

## Load necessary modules 
module purge all
module load perl/5.16

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/bbcar/sample_names_tumor_normal.txt

## Set input and output directories 
din='/projects/b1131/saya/bbcar/02_variant_calls/tumor_normal'
dout='/projects/b1131/saya/bbcar/03_annotated_variants/annovar/tumor_normal'
dav='/projects/b1131/saya/bbcar/tools/annovar'

##########################
#### With generic PON ####
##########################

fin=${input_args[$SLURM_ARRAY_TASK_ID]}_DPfiltered.vcf
fout=${input_args[$SLURM_ARRAY_TASK_ID]}

perl ${dav}/table_annovar.pl \
    ${din}/${fin} \
    ${dav}/humandb/ \
    -buildver hg38 \
    -out $dout/$fout \
    -remove \
    -protocol refGene,knownGene,ensGene,avsnp150,dbnsfp35a,dbnsfp31a_interpro,exac03,gnomad211_exome,gnomad211_genome \
    -operation g,g,g,f,f,f,f,f,f \
    -nastring . -vcfinput
#

########################
#### With BBCAR PON ####
########################

fin=${input_args[$SLURM_ARRAY_TASK_ID]}_DPfiltered_bbcarpon.vcf
fout=${input_args[$SLURM_ARRAY_TASK_ID]}_bbcarpon

perl ${dav}/table_annovar.pl \
    ${din}/${fin} \
    ${dav}/humandb/ \
    -buildver hg38 \
    -out $dout/$fout \
    -remove \
    -protocol refGene,knownGene,ensGene,avsnp150,dbnsfp35a,dbnsfp31a_interpro,exac03,gnomad211_exome,gnomad211_genome \
    -operation g,g,g,f,f,f,f,f,f \
    -nastring . -vcfinput
#
