import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

"""
### Note: First run the below on commandline 
module purge all
module load bcftools/1.10.1

#####################
#### Generic PON ####
#####################

## Tumor only - pre-DP filter 
for fn in $(ls /projects/b1131/saya/bbcar/02_variant_calls/tumor_only/*_filtered.vcf); do
    bcftools stats $fn | grep "number of records:" >> /projects/b1131/saya/bbcar/02_variant_calls/tumor_only/varcounts_preDPfilter.txt;
done

## Tumor only - post-DP filter 
for fn in $(ls /projects/b1131/saya/bbcar/02_variant_calls/tumor_only/*_DPfiltered.vcf); do
    bcftools stats $fn | grep "number of records:" >> /projects/b1131/saya/bbcar/02_variant_calls/tumor_only/varcounts_postDPfilter.txt;
done

## Germline only - pre-DP filter 
for fn in $(ls /projects/b1131/saya/bbcar/02_variant_calls/germline_only/*_haplotype.vcf); do
    bcftools stats $fn | grep "number of records:" >> /projects/b1131/saya/bbcar/02_variant_calls/germline_only/varcounts_preDPfilter.txt;
done

## Germline only - post-DP filter 
for fn in $(ls /projects/b1131/saya/bbcar/02_variant_calls/germline_only/*_DPfiltered.vcf); do
    bcftools stats $fn | grep "number of records:" >> /projects/b1131/saya/bbcar/02_variant_calls/germline_only/varcounts_postDPfilter.txt;
done

## Tumor + normal - pre-DP filter 
for fn in $(ls /projects/b1131/saya/bbcar/02_variant_calls/tumor_normal/*_filtered.vcf); do
    bcftools stats $fn | grep "number of records:" >> /projects/b1131/saya/bbcar/02_variant_calls/tumor_normal/varcounts_preDPfilter.txt;
done

## Tumor + normal - post-DP filter 
for fn in $(ls /projects/b1131/saya/bbcar/02_variant_calls/tumor_normal/*_DPfiltered.vcf); do
    bcftools stats $fn | grep "number of records:" >> /projects/b1131/saya/bbcar/02_variant_calls/tumor_normal/varcounts_postDPfilter.txt;
done

###################
#### BBCAR PON ####
###################

## Tumor only - pre-DP filter 
for fn in $(ls /projects/b1131/saya/bbcar/02_variant_calls/tumor_only/*_filtered_bbcarpon.vcf); do
    bcftools stats $fn | grep "number of records:" >> /projects/b1131/saya/bbcar/02_variant_calls/tumor_only/varcounts_preDPfilter_bbcarpon.txt;
done

## Tumor only - post-DP filter 
for fn in $(ls /projects/b1131/saya/bbcar/02_variant_calls/tumor_only/*_DPfiltered_bbcarpon.vcf); do
    bcftools stats $fn | grep "number of records:" >> /projects/b1131/saya/bbcar/02_variant_calls/tumor_only/varcounts_postDPfilter_bbcarpon.txt;
done

## Tumor + normal - pre-DP filter 
for fn in $(ls /projects/b1131/saya/bbcar/02_variant_calls/tumor_normal/*_filtered_bbcarpon.vcf); do
    bcftools stats $fn | grep "number of records:" >> /projects/b1131/saya/bbcar/02_variant_calls/tumor_normal/varcounts_preDPfilter_bbcarpon.txt;
done

## Tumor + normal - post-DP filter 
for fn in $(ls /projects/b1131/saya/bbcar/02_variant_calls/tumor_normal/*_DPfiltered_bbcarpon.vcf); do
    bcftools stats $fn | grep "number of records:" >> /projects/b1131/saya/bbcar/02_variant_calls/tumor_normal/varcounts_postDPfilter_bbcarpon.txt;
done

"""

######################################
#### Tumor only samples' variants ####
######################################

### Generic PON ###

## Pre-filter 
with open('/projects/b1131/saya/bbcar/02_variant_calls/tumor_only/varcounts_preDPfilter.txt', 'r') as f:
    lines = f.readlines()

pre_to = [] # pre-filter, tumor only 
for line in lines:
    pre_to.append(int(line.strip().split()[-1]))

## Post-filter 
with open('/projects/b1131/saya/bbcar/02_variant_calls/tumor_only/varcounts_postDPfilter.txt', 'r') as f:
    lines = f.readlines()

post_to = [] # pre-filter, tumor only 
for line in lines:
    post_to.append(int(line.strip().split()[-1]))

### BBCAR PON ###

## Pre-filter 
with open('/projects/b1131/saya/bbcar/02_variant_calls/tumor_only/varcounts_preDPfilter_bbcarpon.txt', 'r') as f:
    lines = f.readlines()

pre_to_bbpon = [] # pre-filter, tumor only 
for line in lines:
    pre_to_bbpon.append(int(line.strip().split()[-1]))

## Post-filter 
with open('/projects/b1131/saya/bbcar/02_variant_calls/tumor_only/varcounts_postDPfilter_bbcarpon.txt', 'r') as f:
    lines = f.readlines()

post_to_bbpon = [] # pre-filter, tumor only 
for line in lines:
    post_to_bbpon.append(int(line.strip().split()[-1]))


#########################################
#### Germline only samples' variants ####
#########################################

## Pre-filter 
with open('/projects/b1131/saya/bbcar/02_variant_calls/germline_only/varcounts_preDPfilter.txt', 'r') as f:
    lines = f.readlines()

pre_go = [] # pre-filter, germline only 
for line in lines:
    pre_go.append(int(line.strip().split()[-1]))

## Post-filter 
with open('/projects/b1131/saya/bbcar/02_variant_calls/germline_only/varcounts_postDPfilter.txt', 'r') as f:
    lines = f.readlines()

post_go = [] # pre-filter, germline only 
for line in lines:
    post_go.append(int(line.strip().split()[-1]))


##########################################
#### Tumor + normal samples' variants ####
##########################################

### Generic PON ###

## Pre-filter 
with open('/projects/b1131/saya/bbcar/02_variant_calls/tumor_normal/varcounts_preDPfilter.txt', 'r') as f:
    lines = f.readlines()

pre_tn = [] # pre-filter, tumor normal
for line in lines:
    pre_tn.append(int(line.strip().split()[-1]))

## Post-filter 
with open('/projects/b1131/saya/bbcar/02_variant_calls/tumor_normal/varcounts_postDPfilter.txt', 'r') as f:
    lines = f.readlines()

post_tn = [] # pre-filter, tumor only 
for line in lines:
    post_tn.append(int(line.strip().split()[-1]))

### BBCAR PON ###

## Pre-filter 
with open('/projects/b1131/saya/bbcar/02_variant_calls/tumor_normal/varcounts_preDPfilter_bbcarpon.txt', 'r') as f:
    lines = f.readlines()

pre_tn_bbpon = [] # pre-filter, tumor normal
for line in lines:
    pre_tn_bbpon.append(int(line.strip().split()[-1]))

## Post-filter 
with open('/projects/b1131/saya/bbcar/02_variant_calls/tumor_normal/varcounts_postDPfilter_bbcarpon.txt', 'r') as f:
    lines = f.readlines()

post_tn_bbpon = [] # pre-filter, tumor only 
for line in lines:
    post_tn_bbpon.append(int(line.strip().split()[-1]))


##############
#### Plot ####
##############

plt.violinplot([pre_to, post_to, pre_go, post_go, pre_tn, post_tn], np.arange(6), showextrema=True, showmedians=True, bw_method='silverman')
plt.ylabel('Number of records')
plt.xticks(np.arange(6), ['Tissue-only\n pre-filter', 'Tissue-only\n post-filter', 'Germline-only\n pre-filter', 'Germline-only\n post-filter', 'Tissue-germline\n pre-filter', 'Tissue-germline\n post-filter'], rotation=45, ha='right')
plt.title('Number of variants called per patient - generic PON')
plt.ylim(0, 20000) # 17500
plt.tight_layout()
plt.savefig('/home/srd6051/bbcar_varcounts_DPfilter_violin.png')
plt.close()

plt.violinplot([pre_to_bbpon, post_to_bbpon, pre_go, post_go, pre_tn_bbpon, post_tn_bbpon], np.arange(6), showextrema=True, showmedians=True, bw_method='silverman')
plt.ylabel('Number of records')
plt.xticks(np.arange(6), ['Tissue-only\n pre-filter', 'Tissue-only\n post-filter', 'Germline-only\n pre-filter', 'Germline-only\n post-filter', 'Tissue-germline\n pre-filter', 'Tissue-germline\n post-filter'], rotation=45, ha='right')
plt.title('Number of variants called per patient - BBCAR PON')
plt.ylim(0, 20000) # 17500
plt.tight_layout()
plt.savefig('/home/srd6051/bbcar_varcounts_DPfilter_violin_bbcarpon.png')
plt.close()
