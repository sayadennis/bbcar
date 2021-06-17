import os
import sys
import numpy as np
import pandas as pd

din = '/projects/b1122/saya/04_cleaned_cnv'
dout = '/projects/b1122/saya/scanmap_data'

labdn = '/projects/b1122/saya' # directory where label/target is
ixdn = '/projects/b1122/saya/scanmap_data' # directory where to save train-val-test indices

###############################################
#### First process the gene-level CNV data ####
###############################################

# load data
orig = pd.read_csv(os.path.join(din, 'gene_thres_abs.csv'), header=0, index_col=0)
t1 = pd.read_csv(os.path.join(din, 'gene_thres_thres005_abs.csv'), header=0, index_col=0)
ts1 = pd.read_csv(os.path.join(din, 'gene_thres_thres010_abs.csv'), header=0, index_col=0)

# loop through data and update with new index
for df in [orig, t1, ts1]:
    df['patid'] = None
    for patid in list(df.index):
        int_id = int(patid.split('_')[0])
        df.iloc[[x == patid for x in df.index],-1] = int_id
    df.set_index(['patid'], drop=True, inplace=True)
    df.sort_index(inplace=True)
    df.reset_index(inplace=True, drop=True) # reset in place

# save as pickle
import pickle

outfn = 'gene_thres_abs.pik'
with open(os.path.join(dout, outfn), 'wb') as fh:
    pickle.dump(orig, fh)

outfn = 'gene_thres_thres005_abs.pik'
with open(os.path.join(dout, outfn), 'wb') as fh:
    pickle.dump(orig, fh)

outfn = 'gene_thres_thres010_abs.pik'
with open(os.path.join(dout, outfn), 'wb') as fh:
    pickle.dump(orig, fh)

################################################
#### Next process the region-level CNV data ####
################################################

# load data
orig = pd.read_csv(os.path.join(din, 'reg_thres_conf90.csv'), header=0, index_col=0).T

#### Set row index to be integer values representing patient IDs #### 
orig['patid'] = None
for patid in list(orig.index):
    int_id = int(patid.split('_')[0])
    orig.iloc[[x == patid for x in orig.index],-1] = int_id
orig.set_index(['patid'], drop=True, inplace=True)
orig.sort_index(inplace=True)
orig.reset_index(inplace=True, drop=True) # reset in place

# save as pickle
outfn = 'reg_thres.pik'
with open(os.path.join(dout, outfn), 'wb') as fh:
    pickle.dump(orig, fh)

#################################
#### Also reindex label file ####
#################################

# load data 
labels = pd.read_csv(os.path.join(labdn, 'bbcar_label.csv'), header=0, index_col=0)

# reindex
labels['patid'] = None
for patid in list(labels.index):
    int_id = int(patid.split('_')[0])
    labels.iloc[[x == patid for x in labels.index],-1] = int_id
labels.set_index(['patid'], drop=True, inplace=True)
labels.sort_index(inplace=True)
labels.reset_index(inplace=True, drop=True)

# save as CSV
labels.to_csv(os.path.join(labdn, 'bbcar_label_intindex.csv'), header=True, index=True)

##################################################
#### Store train-val-test indices for ScanMap ####
##################################################

## Store train and test indices for ScanMap 
from sklearn.model_selection import train_test_split
seed = 0
df = orig # the one that has integer index 
X_train, X_test, y_train, y_test = train_test_split(df, labels, test_size=31, random_state=seed, stratify=labels)
X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=30, random_state=seed, stratify=y_train)
train_index = list(X_train.index)
test_index = list(X_test.index)
val_index = list(X_val.index)

# save indices as 
with open(os.path.join(ixdn, 'train_indices_0.15val_0.15te.csv'), 'w') as fh:
    for item in train_index:
        fh.write('%s\n' % item)

with open(os.path.join(ixdn, 'test_indices_0.15val_0.15te.csv'), 'w') as fh:
    for item in test_index:
        fh.write('%s\n' % item)

with open(os.path.join(ixdn, 'val_indices_0.15val_0.15te.csv'), 'w') as fh:
    for item in val_index:
        fh.write('%s\n' % item)

