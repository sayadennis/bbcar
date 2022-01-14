import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

din='/projects/b1122/saya'
dout='/projects/b1122/saya/statistical_tests'

data = pd.read_csv(f'{din}/04_cleaned_cnv/gene_copy_conf90.csv', index_col=0)
labs = pd.read_csv(f'{din}/bbcar_label.csv', index_col=0)

cases = data.iloc[labs.values==1,:]
controls = data.iloc[labs.values==0,:]

amplified = pd.DataFrame(None, columns=['gene', 't stat', 'p-value'])
deleted = pd.DataFrame(None, columns=['gene', 't stat', 'p-value'])

for gene in data.columns:
    t, p = ttest_ind(cases[gene], controls[gene])
    if p <= 0.05:
        if t > 0:
            amplified = pd.concat((amplified, pd.DataFrame.from_dict({'gene' : [gene], 't stat' : [t], 'p-value' : [p]})), ignore_index=True)
        else:
            deleted = pd.concat((deleted, pd.DataFrame.from_dict({'gene' : [gene], 't stat' : [t], 'p-value' : [p]})), ignore_index=True)
            # deleted.append(gene)

amplified.to_csv(f'{dout}/amplified_genes.csv')
deleted.to_csv(f'{dout}/deleted_genes.csv')
