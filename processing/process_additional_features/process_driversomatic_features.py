import os
import sys
import numpy as np
import pandas as pd

din = '/projects/b1122/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Mutect/Somatic_Matrix'
dout = '/projects/b1122/saya/additional_features'

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

sample_ids = list(pd.read_csv('bbcar/all_samples.txt', index_col=0, header=None).index)

# dataframe to record the counts of somatic mutations in each gene
somatic_ct_mx = pd.DataFrame(None, dtype=int, index=sample_ids)

# remove patient ID whose mutational data seem to be missing 
sample_ids_rm467 = sample_ids.copy()
sample_ids_rm467.remove(467) # this patient doesn't seem to have mutational data 

# loop through genes and append number of mutations
for genename in uniq_drivergenes:
    if genename in os.listdir(din):
        counts = pd.read_csv(os.path.join(din, genename, 'number.csv'), index_col=0)
        counts = counts.loc[sample_ids_rm467]
        sum_counts = counts.sum(axis=1)
        # add row for patient 467 and set values to zero
        sum_counts = sum_counts.append(pd.Series(0, index=[467], name=genename))
        sum_counts.sort_index(inplace=True)
        somatic_ct_mx[genename] = sum_counts.values
    else:
        continue


#### Save necessary information in ML-appropriate format #### 
somatic_ct_mx.to_csv(os.path.join(dout, 'bbcar_driversomatic_studyindex.csv'), header=True, index=True)

# save CSV with integer index too
somatic_ct_mx.reset_index(drop=True).to_csv(os.path.join(dout, 'bbcar_driversomatic_intindex.csv'), header=True, index=True)
