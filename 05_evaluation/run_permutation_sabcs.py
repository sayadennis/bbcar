import os
import sys
import numpy as np
import pandas as pd
import glob

sys.path.append('bbcar/src/evaluation')
from RandomPermutation import get_pval

labdn = '/projects/b1122/saya'
ixdir = '/projects/b1122/saya/indices'
resultsdir = 'bbcar/model_performance'

test_indices = list(pd.read_csv(f'{ixdir}/test_indices_0.1val_0.2te.csv', header=None, index_col=0).index)
y_test = pd.read_csv(f'{labdn}/bbcar_label_intindex.csv', header=0, index_col=0).loc[test_indices,:].to_numpy().ravel()

fdict = {
    'clinical' : f'{resultsdir}/noCNV/pred_noCNV_clin_C0.001.csv',
    'mutational signature' : f'{resultsdir}/noCNV/pred_noCNV_mut_C1.csv',
    'factorized CNA' : f'{resultsdir}/results_scanmap_siamese/pred_regthres_nocf_k240_wcls0.01_C0.1_lr0.001.txt',
    'factorized CNA + clinical' : f'{resultsdir}/results_scanmap_siamese/pred_regthres_clin_k120_wcls0.01_C0.01_lr0.01.txt', 
    'factorized CNA + mutational signature' : f'{resultsdir}/results_scanmap_siamese/pred_regthres_mut_k120_wcls0.01_C0.01_lr0.001.txt',
    'factorized CNA + clinical + mutational signature' : f'{resultsdir}/results_scanmap_siamese/pred_regthres_clin_mut_k120_wcls1_C0.01_lr0.001.txt',
    'factorized CNA + clinical + mutational signature + driver somatic mutations' : f'{resultsdir}/results_scanmap_siamese/pred_regthres_clin_mut_driversomatic_k300_wcls1_C0.01_lr0.001.txt',
    'factorized CNA + clinical + mutational signature + PRS' : f'{resultsdir}/results_scanmap_siamese/pred_regthres_clin_mut_prs_k240_wcls1_C0.01_lr0.001.txt'
}

def compare(feat1, feat2):
    print(f'#### {feat1} vs. {feat2} ####')
    feature1 = pd.read_csv(fdict[feat1], header=None, dtype=int)
    feature2 = pd.read_csv(fdict[feat2], header=None, dtype=int)
    pval, pf_feat1, pf_feat2 = get_pval(feature1, feature2, y_test)
    print(f'Performance for {feat1}: {pf_feat1}')
    print(f'Performance for {feat2}: {pf_feat2}')
    print(f'p={pval}\n')


#### Comparisons with clinical as baseline ####

# Clinical vs. CNA
compare('clinical', 'factorized CNA')

# Clinical vs. CNA + clinical 
compare('clinical', 'factorized CNA + clinical')

# Clinical vs. CNA + clinical + mutational signature
compare('clinical', 'factorized CNA + clinical + mutational signature')

#### Comparisons with mutational signature as baseline ####

# Mut sig vs. CNA 
compare('mutational signature', 'factorized CNA')

# Mut sig vs. CNA + mut sig 
compare('mutational signature', 'factorized CNA + mutational signature')

# Mut sig vs. CNA + clinical + mutational signature
compare('mutational signature', 'factorized CNA + clinical + mutational signature')

#### Comparisons with factorized CNA as baseline ####

# CNA vs. CNA + clinical 
compare('factorized CNA', 'factorized CNA + clinical')

# CNA vs. CNA + mut sig 
compare('factorized CNA', 'factorized CNA + mutational signature')

# CNA vs. CNA + clinical + mut sig 
compare('factorized CNA', 'factorized CNA + clinical + mutational signature')

# CNA vs. CNA + clinical + mut sig + PRS 
compare('factorized CNA', 'factorized CNA + clinical + mutational signature + PRS')

# CNA vs. CNA + clinical + mut sig + driver somatic mutations 
compare('factorized CNA', 'factorized CNA + clinical + mutational signature + driver somatic mutations')
