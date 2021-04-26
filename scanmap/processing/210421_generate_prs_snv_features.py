# The purpose of this script is to create additional confounding features for ScanMap 
# which represent the presence of risk SNVs described in a polygenic risk scoring paper (Mavaddat et al. 2019)

import os
import sys
import numpy as np
import pandas as pd

din = '/share/fsmresfiles/bbcar/additional_features'
dout = '/share/fsmresfiles/bbcar/scanmap_data'

prs = pd.read_csv(os.path.join(din, 'bbcar_prs_snp.csv'), index_col=0)
target = pd.read_csv('/share/fsmresfiles/bbcar/bbcar_label_studyidindex.csv', index_col=0)

## Remove duplicate rows and align patient IDs with target 
prs = prs[~prs.index.duplicated(keep='first')]
prs = prs.loc[target.index,:]

## Remove SNVs where all patients are zero or number of patients who have that SNV is <= 3 
rmzero = prs.iloc[:,[x!=0 for x in prs.sum().values]]
rmthree = prs.iloc[:,[x<=3 for x in prs.sum().values]]

## Save features by themselves
rmzero.reset_index(inplace=True)
rmzero.to_csv(os.path.join(dout, 'prs_snp_intindex.csv'), header=True, index=True)

rmthree.reset_index(inplace=True)
rmthree.to_csv(os.path.join(dout, 'prs_snp_rm3_intindex.csv'), header=True, index=True)

## Save features concatenated with clinical or clinical + mutational signatures 
clin = pd.read_csv(os.path.join(dout, 'pt_demo_intindex.csv'), index_col=0)
clinmut = pd.read_csv(os.path.join(dout, 'pt_demo_clinmut_intindex.csv'), index_col=0)

clin_rmzero = pd.concat((clin, rmzero))
clinmut_rmzero = pd.concat((clinmut, rmzero))

clin_rmzero.to_csv(os.path.join(dout, 'clin_prssnp.csv'), header=True, index=True)
clinmut_rmzero.to_csv(os.path.join(dout, 'clinmut_prssnp.csv'), header=True, index=True)

clin_rmthree = pd.concat((clin, rmthree))
clinmut_rmthree = pd.concat((clinmut, rmthree))

clin_rmthree.to_csv(os.path.join(dout, 'clin_prssnp_rm3.csv'), header=True, index=True)
clinmut_rmthree.to_csv(os.path.join(dout, 'clinmut_prssnp_rm3.csv'), header=True, index=True)
