import os
import sys
import numpy as np
import pandas as pd
import getopt
import datetime

from sklearn.metrics import roc_auc_score, accuracy_score, precision_score, recall_score, f1_score, balanced_accuracy_score

sys.path.append('bbcar/src/modeling/exploratory_models')
from BBCarMLP import BBCarMLP

import torch
import torch.nn as nn

import importlib
import warnings
warnings.filterwarnings(action='ignore')

opts, extraparams = getopt.getopt(sys.argv[1:], 'g:l:o:i:n:', 
                                  ['genomic=', 'label=', 'outdir=', 'index=', 'n_iter='])

print(sys.argv)
devstr = 'cuda'
config = 'please specify your own configuration string describing, e.g., germline mutations, filtering thresholds in pre-processing steps'
niter = 2000
# lr = 0.01
seed = 1

for o,p in opts:
    if o in ['-g', '--genomic']:
        genfn = p
    if o in ['-l', '--label']:
        labfn = p
    if o in ['-o', '--outdir']:
        outdir = p
    if o in ['-i', '--index']:
        indexdir = p
    if o in ['-n', '--n_iter']:
        niter = int(p)


genomic = pd.read_csv(genfn, index_col=0) # '/projects/b1122/saya/06_modified_data/reg_copy_conf90_intindex.csv'
y = pd.read_csv(labfn, header=0, index_col=0) # '/projects/b1122/saya/bbcar_label_intindex.csv'

if genfn.startswith('reg_copy'):
    feature_type = 'regcopy'
elif genfn.startswith('reg_thres'):
    feature_type = 'regthres'
elif genfn.startswith('gene_copy'):
    feature_type = 'genecopy'
elif genfn.startswith('gene_thres'):
    feature_type = 'genethres'
else:
    print('Unknown feature type name for file %s' % genfn)

# indexdir = '/projects/b1122/saya/indices'
train_indices = pd.read_csv('%s/train_indices_0.1val_0.2te.csv' % (indexdir), header=None) #  _5run
test_indices = pd.read_csv('%s/test_indices_0.1val_0.2te.csv' % (indexdir), header=None)
val_indices = pd.read_csv('%s/val_indices_0.1val_0.2te.csv' % (indexdir), header=None)

X = np.array(genomic)
y, yuniques = pd.factorize(y.values.ravel(), sort=True)

r = 0
X = X.astype(float)
device = torch.device(devstr)
devstr = 'cuda'

train_index = train_indices[r]; val_index = val_indices[r]; test_index = test_indices[r]
X_train, X_val, X_test = X[train_index], X[val_index], X[test_index]
y_train, y_val, y_test = y[train_index], y[val_index], y[test_index]

print('C,lr,n_iter,tr bal acc,val bal acc,te acc,te bal acc,te precis,te recall,te f1,celoss') # ,best iter,mse,mse tr,mse val,mse te
for lr in [0.0001, 0.001, 0.01]: # 0.00001, 0.0001 
    for C in [0.001, 0.01, 0.1, 1, 10, 100, 1000]: #1000, 0.001 
        fn = '%s/bbcarmlp_%s_C%s_lr%s_niter%s_seed%s.p' % (outdir, feature_type, C, lr, niter, seed) # scanmap%d/s%d/ # niter, seed, 
        m = BBCarMLP(X_train, X_val, y_train, y_val, n_iter=niter, fn=fn, C=C, lr=lr, device=device)
        m.fit()

        chkpt = torch.load(fn)
        m.load_state_dict(chkpt['state_dict'])
        accval = chkpt['best_val_acc']

        m.eval()

        y_tr_pred = m.predict(X_train)
        y_te_pred = m.predict(X_test)
        acctr = balanced_accuracy_score(y_train, y_tr_pred)
        accte = accuracy_score(y_test, y_te_pred)
        balaccte = balanced_accuracy_score(y_test, y_te_pred)
        preciste = precision_score(y_test, y_te_pred)
        recallte = recall_score(y_test, y_te_pred)
        f1te = f1_score(y_test, y_te_pred)
        vce = chkpt['celoss']

        print('%s,%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f' % 
            (C, lr, niter, acctr, accval, accte, balaccte, preciste, recallte, f1te, vce))
