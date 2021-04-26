# The purpose of this script is to create additional confounding features for ScanMap 
# which represent the number of driver mutations in genes known to be mutated in samples with similar mutational signatures as our BBCar samples

import os
import sys
import numpy as np
import pandas as pd

din = '/projects/b1122/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/Somatic_Matrix'
dout = '/projects/b1122/saya/bbcar_project/additional_features'

driver_genes = [
    'CTNNB1', 'LINC00290', 'ALB', # SigD genes
    'BRCA2', 'TP53', 'BRCA1', 'MYC', 'ARID1A', 'LINC00290', 'NF1', # SigK genes
    'TP53', # SigG genes
    'CDKN2A', 'CDKN2B', 'TP53', 'KRAS', 'EGFR', 'SMAD4', 'APC', 'BRD4'
]

# remove duplicates
uniq_drivergenes = []
for gene in driver_genes:
    if gene not in uniq_drivergenes:
        uniq_drivergenes.append(gene)

# list of patient IDs (sorted as integers)
with open('bbcar_project/patient_ids_with_cnv.txt', 'r') as f:
    lines = f.readlines()

patids = []
for line in lines:
    patids.append(int(line.rstrip()))

# dataframe to record the counts of somatic mutations in each gene
somatic_ct_mx = pd.DataFrame(None, dtype=int, index=patids)

# remove patient ID whose mutational data seem to be missing 
patids_rm467 = patids.copy()
patids_rm467.remove(467) # this patient doesn't seem to have mutational data 

# loop through genes and append number of mutations
for genename in uniq_drivergenes:
    if genename in os.listdir(din):
        counts = pd.read_csv(os.path.join(din, genename, 'number.csv'), index_col=0)
        counts = counts.loc[patids_rm467]
        sum_counts = counts.sum(axis=1)
        # add row for patient 467 and set values to zero
        sum_counts = sum_counts.append(pd.Series(0, index=[467], name=genename))
        sum_counts.sort_index(inplace=True)
        somatic_ct_mx[genename] = sum_counts.values
    else:
        continue

somatic_ct_mx.to_csv(os.path.join(dout, 'drivergenes_somatic_mut_ct.csv'))

somatic_ct_mx.reset_index(inplace=True, drop=True)
somatic_ct_mx.to_csv(os.path.join(dout, 'drivergenes_somatic_mut_ct_intindex.csv'))
