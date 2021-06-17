import os
import sys
import pickle
import numpy as np
import pandas as pd

sys.path.append('bbcar/src')
from applyThres import applyThres

genfn = '/projects/b1122/saya/scanmap_data/gene_thres_abs.pik'

f = open(genfn, 'rb')
genomic = pickle.load(f)
f.close()

red005 = applyThres(genomic, thres=0.05, twotail=True)
print('Shape of threshold-applied (a=0.05) matrix: {}x{}'.format(red005.shape[0], red005.shape[1]))

red010 = applyThres(genomic, thres=0.10, twotail=True)
print('Shape of threshold-applied (a=0.10) matrix: {}x{}'.format(red010.shape[0], red010.shape[1]))

with open('/projects/b1122/saya/scanmap_data/gene_thres_thres005_abs.pik', 'wb') as f:
    pickle.dump(red005, f)

with open('/projects/b1122/saya/scanmap_data/gene_thres_thres010_abs.pik', 'wb') as f:
    pickle.dump(red010, f)
