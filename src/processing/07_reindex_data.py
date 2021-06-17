import os
import sys
import numpy as np
import pandas as pd

din = '/projects/b1122/saya/04_cleaned_cnv'
dout = '/projects/b1122/saya/06_modified_data'

labdn = '/projects/b1122/saya' # directory where label/target is
ixdn = '/projects/b1122/saya/indices' # directory where to save train-val-test indices

#############################################################
#### First process the original, cleaned GISTIC2 outputs ####
#############################################################

datadict = {}
flist = [
    'gene_copy_conf90.csv',
    'gene_thres_conf90.csv',
    'reg_copy_conf90.csv',
    'reg_thres_conf90.csv'
]

# Read the four files into a dictionary { 'filename' : pd.DataFrame }
for fn in flist:
    if 'gene_' in fn:
        datadict[fn] = pd.read_csv(os.path.join(din, fn), index_col=0)
    elif 'reg_' in fn:
        datadict[fn] = pd.read_csv(os.path.join(din, fn), index_col=0).T
    else:
        print('Unknown pattern: %s' % fn)

# loop through data and update with new index
for fn in datadict.keys():
    datadict[fn]['patid'] = None
    for patid in list(datadict[fn].index):
        int_id = int(patid.split('_')[0])
        datadict[fn].iloc[[x == patid for x in datadict[fn].index],-1] = int_id
    datadict[fn].set_index(['patid'], drop=True, inplace=True)
    datadict[fn].sort_index(inplace=True)
    datadict[fn].to_csv(os.path.join(dout, fn.split('.')[0] + '_studyindex.csv'))
    datadict[fn].reset_index(inplace=True, drop=True) # reset in place
    datadict[fn].to_csv(os.path.join(dout, fn.split('.')[0] + '_intindex.csv'))


##########################################################################
#### Next process the frequency-based feature-reduced gene-thres data ####
##########################################################################

datadict = {}
flist = [
    'gene_thres_abs.csv',
    'gene_thres_thres005_abs.csv',
    'gene_thres_thres010_abs.csv'
]

# Read the three files into a dictionary { 'filename' : pd.DataFrame }
for fn in flist:
    datadict[fn] = pd.read_csv(os.path.join(dout, fn), header=0, index_col=0)

# loop through data and update with new index
for fn in datadict.keys():
    datadict[fn]['patid'] = None
    for patid in list(datadict[fn].index):
        int_id = int(patid.split('_')[0])
        datadict[fn].iloc[[x == patid for x in datadict[fn].index],-1] = int_id
    datadict[fn].set_index(['patid'], drop=True, inplace=True)
    datadict[fn].sort_index(inplace=True)
    datadict[fn].to_csv(os.path.join(dout, fn.split('.')[0] + '_studyindex.csv'))
    datadict[fn].reset_index(inplace=True, drop=True) # reset in place
    datadict[fn].to_csv(os.path.join(dout, fn.split('.')[0] + '_intindex.csv'))

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

# save with study ID as index
labels.to_csv(os.path.join(labdn, 'bbcar_label_studyindex.csv'), header=True, index=True)

# save with sequential integers as index
labels.reset_index(inplace=True, drop=True)
labels.to_csv(os.path.join(labdn, 'bbcar_label_intindex.csv'), header=True, index=True)

##################################################
#### Store train-val-test indices for ScanMap ####
##################################################

## Store train and test indices for ScanMap 
from sklearn.model_selection import train_test_split
seed = 0
df = datadict['gene_thres_abs.csv'] # now has integer index 
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

