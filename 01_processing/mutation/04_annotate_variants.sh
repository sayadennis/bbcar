#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH --array=0-239
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 24:00:00
#SBATCH --mem=8G
#SBATCH --mail-user=sayarenedennis@northwestern.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=annotate
#SBATCH --output=/projects/b1131/saya/new_bbcar/out/annotate_variants_%a.out

cd /projects/b1131/saya/new_bbcar/

## Load necessary modules 
module purge all
module load perl/5.16
module load java/jdk11.0.10

## Define input arguments for job array 
IFS=$'\n' read -d '' -r -a input_args < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_tissue.txt
sampleid=${input_args[$SLURM_ARRAY_TASK_ID]}

## Obtain IDs of samples that have germline
IFS=$'\n' read -d '' -r -a tissue < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_tissue.txt
IFS=$'\n' read -d '' -r -a germline < /projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_germline.txt

data_dir="/projects/b1131/saya/new_bbcar/data/02a_mutation"
dav='/projects/b1131/saya/bbcar/tools/annovar'

function run_annovar() {
    local sampleid=$1
    local mode=$2  # tissue_only etc.

    if [ "$mode" = "germline_only" ]; then
        ext="DPfiltered"
    else
        ext="DPfiltered_classicalAF"
    fi

    ## Set input and output directories 
    din="${data_dir}/02_variant_calls/${mode}"
    dout="${data_dir}/03_annotated_variants/annovar/${mode}"

    mkdir -p $dout
    
    fin=${sampleid}_${ext}.vcf
    fout=${sampleid}

    perl ${dav}/table_annovar.pl \
        ${din}/${fin} \
        ${dav}/humandb/ \
        -buildver hg38 \
        -out $dout/$fout \
        -remove \
        -protocol refGene,knownGene,ensGene,avsnp150,dbnsfp35a,dbnsfp31a_interpro,exac03,gnomad211_exome \
        -operation g,g,g,f,f,f,f,f \
        -nastring . -vcfinput
}

function run_snpeff() {
    local sampleid=$1
    local mode=$2

    if [ "$mode" = "germline_only" ]; then
        ext="DPfiltered"
    else
        ext="DPfiltered_classicalAF"
    fi

    ## Set input and output directories 
    din="${data_dir}/02_variant_calls/${mode}"
    dout="${data_dir}/03_annotated_variants/snpeff/${mode}"
    dspf='/projects/b1131/saya/bbcar/tools/snpEff'

    mkdir -p $dout
    
    java -Xmx5120m -jar ${dspf}/snpEff.jar hg38 \
        ${din}/${sampleid}_${ext}.vcf \
        -stats ${dout}/${sampleid}_snpEff_genes.txt > \
        ${dout}/${sampleid}_snpeff.vcf
}

function run_vep() {
    local sampleid=$1
    local mode=$2

    if [ "$mode" = "germline_only" ]; then
        ext="DPfiltered"
    else
        ext="DPfiltered_classicalAF"
    fi

    ## Set input and output directories 
    din="${data_dir}/02_variant_calls/${mode}"
    dout="${data_dir}/03_annotated_variants/vep/${mode}"
    dvep='/projects/b1131/saya/bbcar/tools/ensembl-vep'
    
    mkdir -p $dout
    
    ${dvep}/vep \
        -i ${din}/${sampleid}_${ext}.vcf \
        -o ${dout}/${sampleid}_vep.vcf \
        --force_overwrite \
        -offline --cache --dir /projects/b1131/saya/bbcar/tools/.vep
}

## Tissue only
if [[ " ${tissue[*]} " =~ " ${sampleid} " ]]; then
    echo "########## Annotating variants on tissue sample ${sampleid} ##########"
    run_annovar ${sampleid} 'tissue_only'
    run_snpeff ${sampleid} 'tissue_only'
    run_vep ${sampleid} 'tissue_only'
else
    echo '########## Patient ID ${sampleid} is not in the tissue sample list ##########'
fi

## Tissue-normal and germline-only
if [[ " ${germline[*]} " =~ " ${sampleid} " ]]; then
    echo "########## Annotating variants on tissue-germline pair of sample ${sampleid} ##########"
    run_annovar ${sampleid} 'tissue_normal'
    run_snpeff ${sampleid} 'tissue_normal'
    run_vep ${sampleid} 'tissue_normal'
    echo "########## Annotating variants on germline sample ${sampleid} ##########"
    run_annovar ${sampleid} 'germline_only'
    run_snpeff ${sampleid} 'germline_only'
    run_vep ${sampleid} 'germline_only'
else
    echo '########## Patient ID ${sampleid} is not in the germline sample list ##########'
fi

