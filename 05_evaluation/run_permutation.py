import os
import sys
import numpy as np
import pandas as pd
import glob

sys.path.append('bbcar/src/evaluation')
from RandomPermutation import get_pval

datadir = '/projects/b1122/saya/scanmap_data'
ixdir = '/projects/b1122/saya/scanmap_data'
resultsdir = '/home/srd6051/bbcar/scanmap/results'

test_indices = list(pd.read_csv(os.path.join(ixdir, 'test_indices_0.1val_0.2te.csv'), header=None, index_col=0).index)
y_test = pd.read_csv(os.path.join(datadir, 'bbcar_label_intindex.csv'), header=0, index_col=0).loc[test_indices,:].to_numpy().ravel()

cfs = [
    'clin',
    'clin_driversomatic',
    'clin_mut',
    'clin_prs',
    'clin_mut_driversomatic',
    'clin_mut_prs',
    'clin_prs_driversomatic'
]

for cf in cfs:
    flist = glob.glob(resultsdir + '/pred_*_%s_[k|C]*.csv' % cf) # 3 fn with same confounding features and one of {noCNV, genethres, regthres}
    if len(flist) == 2:
        print('Could not find three files for confounding feature: %s. Proceeding with only gene-level CNV features\n' % cf)
        gene_fn = flist[np.where(['genethres' in x for x in flist])[0][0]]
        nocnv_fn = flist[np.where(['noCNV' in x for x in flist])[0][0]]
        genethres = pd.read_csv(gene_fn, header=None, dtype=int)
        nocnv = pd.read_csv(nocnv_fn, header=None, dtype=int)
        print('%s vs. %s' % (gene_fn, nocnv_fn))
        pval, _ = get_pval(genethres, nocnv, y_test)
        print('p=%s' % pval)
        print('\n')
    else:
        gene_fn = flist[np.where(['genethres' in x for x in flist])[0][0]]
        reg_fn = flist[np.where(['regthres' in x for x in flist])[0][0]]
        nocnv_fn = flist[np.where(['noCNV' in x for x in flist])[0][0]]
        genethres = pd.read_csv(gene_fn, header=None, dtype=int)
        regthres = pd.read_csv(reg_fn, header=None, dtype=int)
        nocnv = pd.read_csv(nocnv_fn, header=None, dtype=int)
        print('%s vs. %s' % (gene_fn, nocnv_fn))
        pval, _, _ = get_pval(genethres, nocnv, y_test)
        print('p=%s' % pval)
        print('')
        print('%s vs. %s' % (reg_fn, nocnv_fn))
        pval, _, _ = get_pval(regthres, nocnv, y_test)
        print('p=%s' % pval)
        print('\n')
