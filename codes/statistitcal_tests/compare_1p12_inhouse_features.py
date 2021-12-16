import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

din='/projects/b1122/saya/inhouse_cnv_features'
# dout='/projects/b1122/saya/statistical_tests'

data = pd.read_csv(f'{din}/cytoband_matrix.csv', index_col=0).sort_index()
labs = pd.read_csv(f'{din}/bbcar_label.csv', index_col=0).sort_index()

cases = data.iloc[labs.values==1,:]
controls = data.iloc[labs.values==0,:]

cyto = '1p12'
t, p = ttest_ind(cases[cyto], controls[cyto])
print('#### Results for cytoband 1p12 ####')
print(f't={t}, p={p:.3f}')
print('')

#### gene-level data ####

data = pd.read_csv(f'{din}/refGene_cnv_gene_matrix.csv', index_col=0).sort_index()

cases = data.iloc[labs.values==1,:]
controls = data.iloc[labs.values==0,:]

gene = 'PHGDH'
t, p = ttest_ind(cases[gene], controls[gene])
print(f'#### Results for gene {gene} ####')
print(f't={t}, p={p:.3f}')
print('')

#### gene-level data with GISTIC features ####

data = pd.read_csv('/projects/b1122/saya/06_modified_data/all/gene_copy_conf90_all_studyindex.csv', index_col=0).sort_index()
labs = pd.read_csv('/projects/b1122/saya/bbcar_label_studyid.csv', index_col=0).sort_index()

if not np.all(data.index==labs.index):
    print('!!!! INDEX BETWEEN DATA AND LABELS DO NOT MATCH !!!!')

print('\n...Results with GISTIC features...')
cases = data.iloc[labs.values==1,:]
controls = data.iloc[labs.values==0,:]

gene = 'PHGDH'
t, p = ttest_ind(cases[gene], controls[gene])
print(f'#### Results for gene {gene} ####')
print(f't={t}, p={p:.3f}')
print('')

