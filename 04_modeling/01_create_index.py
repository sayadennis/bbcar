import os
import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

din = '/projects/b1122/saya/01_gatk_analyzed_segments'
dout = '/projects/b1122/saya/indices'

flist = os.listdir(din) # element example: '1004_Tissue.csv' 
studynames = [x.split('.')[0] for x in flist] # element example: '1004_Tissue'
outcome = list(1 if (y=='Tissue') else 0 for y in [x.split('_')[1] for x in studynames])
studyids = [int(x.split('_')[0]) for x in studynames] # element example: 1004 

ix_table = pd.DataFrame([flist, studynames, studyids, outcome], index=['filename', 'samplename', 'studyid', 'outcome']).T
ix_table = ix_table.sort_values('studyid').reset_index(drop=True) # sort by study ID, reset index

ix_table.to_csv(f'{dout}/index_dictionary.csv')

###############################################################
#### Create train-test indices for cross-validation tuning ####
###############################################################

seed = 21
# First, just a train-test split for cross-validation analyses 
train, test = train_test_split(ix_table, test_size=0.2, random_state=seed, shuffle=True, stratify=ix_table['outcome'].values)

# save study IDs
with open(f'{dout}/train_studyid.csv', 'w') as fh:
    for item in list(train['studyid'].values):
        fh.write('%s\n' % item)

with open(f'{dout}/test_studyid.csv', 'w') as fh:
    for item in list(test['studyid'].values):
        fh.write('%s\n' % item)

# save indices
with open(f'{dout}/train_ix.csv', 'w') as fh:
    for item in list(train.index):
        fh.write('%s\n' % item)

with open(f'{dout}/test_ix.csv', 'w') as fh:
    for item in list(test.index):
        fh.write('%s\n' % item)

# Next, another validation split for train-test-val split analyses 
train, test = train_test_split(ix_table, test_size=31, random_state=seed, shuffle=True, stratify=ix_table['outcome'].values)
train, val = train_test_split(train, test_size=30, random_state=seed, shuffle=True, stratify=train['outcome'].values)

# save indices 
with open(f'{dout}/train_indices_0.15val_0.15te.csv', 'w') as fh:
    for item in list(train.index):
        fh.write('%s\n' % item)

with open(f'{dout}/test_indices_0.15val_0.15te.csv', 'w') as fh:
    for item in list(test.index):
        fh.write('%s\n' % item)

with open(f'{dout}/val_indices_0.15val_0.15te.csv', 'w') as fh:
    for item in list(val.index):
        fh.write('%s\n' % item)
