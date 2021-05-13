import os
import sys
import numpy as np
import pandas as pd

dn = '/projects/b1122/saya/additional_features'

clin = pd.read_csv(os.path.join(dn, 'bbcar_clinical_intindex.csv'), index_col=0)
mut = pd.read_csv(os.path.join(dn, 'bbcar_mutsig_intindex.csv'), index_col=0)
driv = pd.read_csv(os.path.join(dn, 'bbcar_driversomatic_intindex.csv'), index_col=0)
prs = pd.read_csv(os.path.join(dn, 'bbcar_prs_intindex.csv'), index_col=0)

## Clinical + mutational signature ## 
pd.concat((clin, mut), axis=1).to_csv(os.path.join(dn, 'bbcar_clin_mut_intindex.csv'))

## Clinical + driver somatic mutations ## 
pd.concat((clin, driv), axis=1).to_csv(os.path.join(dn, 'bbcar_clin_driversomatic_intindex.csv'))

## Clinical + PRS ## 
pd.concat((clin, prs), axis=1).to_csv(os.path.join(dn, 'bbcar_clin_prs_intindex.csv'))

## Clinical + mutational signature + PRS ## 
pd.concat((clin, mut, prs), axis=1).to_csv(os.path.join(dn, 'bbcar_clin_mut_prs_intindex.csv'))

## Clinical + mutational signature + driver somatic mutations ## 
pd.concat((clin, mut, driv), axis=1).to_csv(os.path.join(dn, 'bbcar_clin_mut_driversomatic_intindex.csv'))

## Clinical + PRS + driver somatic mutations ## 
pd.concat((clin, prs, driv), axis=1).to_csv(os.path.join(dn, 'bbcar_clin_prs_driversomatic_intindex.csv'))
