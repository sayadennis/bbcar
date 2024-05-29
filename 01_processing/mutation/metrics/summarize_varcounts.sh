#!/bin/bash

### Note: First run the below on commandline 
# module purge all
module load bcftools/1.10.1

var_dir="/projects/b1131/saya/new_bbcar/data/02a_mutation/02_variant_calls"

#####################
#### Tissue only ####
#####################

## Tissue only - pre-DP filter 
rm ${var_dir}/tissue_only/varcounts_preDPfilter.txt
for fn in $(ls ${var_dir}/tissue_only/*_filtered.vcf); do
    bcftools stats $fn | grep "number of records:" >> ${var_dir}/tissue_only/varcounts_preDPfilter.txt;
done

## Tissue only - post-DP filter 
rm ${var_dir}/tissue_only/varcounts_postDPfilter.txt
for fn in $(ls ${var_dir}/tissue_only/*_DPfiltered.vcf); do
    bcftools stats $fn | grep "number of records:" >> ${var_dir}/tissue_only/varcounts_postDPfilter.txt;
done

## Tissue only - post-FFPE filter 
rm ${var_dir}/tissue_only/varcounts_postFFPEfilter.txt
for fn in $(ls ${var_dir}/tissue_only/*_DPfiltered_classicalAF.vcf); do
    bcftools stats $fn | grep "number of records:" >> ${var_dir}/tissue_only/varcounts_postFFPEfilter.txt;
done

## Tissue only - predicted somatic
rm ${var_dir}/tissue_only/varcounts_predicted_somatic.txt
for fn in $(ls /projects/b1131/saya/new_bbcar/data/02a_mutation/07_predicted_somatic/vcfs/*.vcf); do
    bcftools stats $fn | grep "number of records:" >> ${var_dir}/tissue_only/varcounts_predicted_somatic.txt;
done

#######################
#### Tissue normal ####
#######################

## Tissue + normal - pre-DP filter 
rm ${var_dir}/tissue_normal/varcounts_preDPfilter.txt
for fn in $(ls ${var_dir}/tissue_normal/*_filtered.vcf); do
    bcftools stats $fn | grep "number of records:" >> ${var_dir}/tissue_normal/varcounts_preDPfilter.txt;
done

## Tissue + normal - post-DP filter 
rm ${var_dir}/tissue_normal/varcounts_postDPfilter.txt
for fn in $(ls ${var_dir}/tissue_normal/*_DPfiltered.vcf); do
    bcftools stats $fn | grep "number of records:" >> ${var_dir}/tissue_normal/varcounts_postDPfilter.txt;
done

## Tissue + normal - post-FFPE filter 
rm ${var_dir}/tissue_normal/varcounts_postFFPEfilter.txt
for fn in $(ls ${var_dir}/tissue_normal/*_DPfiltered_classicalAF.vcf); do
    bcftools stats $fn | grep "number of records:" >> ${var_dir}/tissue_normal/varcounts_postFFPEfilter.txt;
done

#######################
#### Germline only ####
#######################

## Germline - pre-DP filter 
rm ${var_dir}/germline_only/varcounts_preDPfilter.txt
for fn in $(ls ${var_dir}/germline_only/*_haplotype.vcf); do
    bcftools stats $fn | grep "number of records:" >> ${var_dir}/germline_only/varcounts_preDPfilter.txt;
done

## Germline - post-DP filter 
rm ${var_dir}/germline_only/varcounts_postDPfilter.txt
for fn in $(ls ${var_dir}/germline_only/*_DPfiltered.vcf); do
    bcftools stats $fn | grep "number of records:" >> ${var_dir}/germline_only/varcounts_postDPfilter.txt;
done

