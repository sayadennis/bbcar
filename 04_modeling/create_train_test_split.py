import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

din = '/projects/b1131/saya/bbcar/data/clinical'
dout = '/projects/b1131/saya/bbcar/train_test_splits'

target = pd.read_csv(f'{din}/bbcar_redcap_label_studyid.csv', index_col=0)

###############################################################
#### Create train-test indices for cross-validation tuning ####
###############################################################

seed = 21
# First, just a train-test split for cross-validation analyses 
train, test = train_test_split(target, test_size=0.2, random_state=seed, shuffle=True, stratify=target.values)

# save study IDs
with open(f'{dout}/train_ix.csv', 'w') as fh:
    for item in list(train.index):
        fh.write('%s\n' % item)

with open(f'{dout}/test_ix.csv', 'w') as fh:
    for item in list(test.index):
        fh.write('%s\n' % item)

# # save indices
# with open(f'{dout}/train_ix.csv', 'w') as fh:
#     for item in list(train.index):
#         fh.write('%s\n' % item)

# with open(f'{dout}/test_ix.csv', 'w') as fh:
#     for item in list(test.index):
#         fh.write('%s\n' % item)

# # Next, another validation split for train-test-val split analyses 
# train, test = train_test_split(ix_table, test_size=31, random_state=seed, shuffle=True, stratify=ix_table['outcome'].values)
# train, val = train_test_split(train, test_size=30, random_state=seed, shuffle=True, stratify=train['outcome'].values)

# # save indices 
# with open(f'{dout}/train_indices_0.15val_0.15te.csv', 'w') as fh:
#     for item in list(train.index):
#         fh.write('%s\n' % item)

# with open(f'{dout}/test_indices_0.15val_0.15te.csv', 'w') as fh:
#     for item in list(test.index):
#         fh.write('%s\n' % item)

# with open(f'{dout}/val_indices_0.15val_0.15te.csv', 'w') as fh:
#     for item in list(val.index):
#         fh.write('%s\n' % item)
