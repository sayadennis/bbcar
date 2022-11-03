import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

datadir = '/projects/b1131/saya/bbcar/data'

cts = pd.read_csv(f'{datadir}/02a_mutation/08_feature_matrix/cts_per_gene.csv', index_col=0)
lab = pd.read_csv(f'{datadir}/clinical/bbcar_redcap_label_studyid.csv', index_col=0)

cosmic_top20 = [
    'PIK3CA',
    'TP53',
    'CDH1',
    'MED12',
    'ESR1',
    'GATA3',
    'KMT2C',
    'MAP3K1',
    'PTEN',
    'LRP1B',
    'ERBB4',
    'ZFHX3',
    'NF1',
    'ERBB2',
    'AKT1',
    'ARID1A',
    'PTPRT',
    'ALK',
    'RUNX1',
    'GRIN2A',
]

ct = 0
for item in cosmic_top20:
    if item not in cts.columns:
        ct += 1

