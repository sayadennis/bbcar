import numpy as np
import pandas as pd

###################################################
#### Create (sample) x (variant) binary matrix ####
###################################################

# Define directories 
dfeat = '/projects/b1131/saya/bbcar/data/02a_mutation/04_ml_features'
dpred = '/projects/b1131/saya/bbcar/data/02a_mutation/07_predicted_somatic'
dout = '/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix'

# Load information on which sample has which variant 
all_var = pd.read_csv(f'{dfeat}/04_imputed/features_imputed.csv')

# Load information on whether mutation is somatic 
y_matched = pd.read_csv(f'{dfeat}/target_matched.csv', index_col=0)
y_nonmatched = pd.read_csv(f'{dpred}/nonmatched.csv', index_col=0)

# Create empty matrix 
X = pd.DataFrame(None, index=all_var.sample_id.unique(), columns=all_var.var_id.unique())

# First, populate the paired calls ("ground truth" from tumor_normal and germline_only)
soma = all_var.iloc[all_var.source.values=='tumor_normal'].groupby(['sample_id', 'var_id']).size()
germ = all_var.iloc[all_var.source.values=='germline_only'].groupby(['sample_id', 'var_id']).size()

for sample_id in soma.index.get_level_values('sample_id'):
    X.loc[sample_id, soma.loc[sample_id].index] = 1

for sample_id in germ.index.get_level_values('sample_id'):
    X.loc[sample_id, germ.loc[sample_id].index] = 0

# Next, populate the predictions - only where entry is NaN
unk = all_var.iloc[all_var.source.values=='tumor_only']
pred_soma = unk.iloc[[bool(y_nonmatched.loc[var_id,'somatic']) for var_id in unk.var_id.values],:]
pred_soma = pred_soma.groupby(['sample_id', 'var_id']).size()

for sample_id in pred_soma.index.get_level_values('sample_id'):
    X.loc[sample_id, pred_soma.loc[sample_id].index] = 1

# Finally, fill NaN with 0 
X = X.fillna(0)

X.to_csv(f'{dout}/binary.csv')
