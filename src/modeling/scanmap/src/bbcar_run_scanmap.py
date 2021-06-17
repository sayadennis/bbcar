import os
import sys
import numpy as np
import pandas as pd
import pickle
import getopt
import datetime
from math import floor, ceil
from scipy.sparse import coo_matrix

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, accuracy_score, precision_score, recall_score, f1_score, balanced_accuracy_score
sys.path.append('bbcar/scanmap/src')
from ScanMap import ScanMap
import torch
import torch.nn as nn

import importlib
import warnings
warnings.filterwarnings(action='ignore')

# print("Starting run ... {}".format(datetime.datetime.now()))

opts, extraparams = getopt.getopt(sys.argv[1:], 'g:c:l:o:i:n:', 
                                  ['genomic=', 'confound=', 'label=', 'outdir=', 'index=', 'n_iter='])

print(sys.argv)
devstr = 'cuda'
config = 'please specify your own configuration string describing, e.g., germline mutations, filtering thresholds in pre-processing steps'
niter = 2000
# lr = 0.01
seed = 1

for o,p in opts:
    if o in ['-g', '--genomic']:
        genfn = p
    if o in ['-c', '--confound']:
        ptsfn = p
    if o in ['-l', '--label']:
        labfn = p
    if o in ['-o', '--outdir']:
        outdir = p
    if o in ['-i', '--index']:
        indexdir = p
    if o in ['-n', '--n_iter']:
        niter = int(p)

f = open(genfn, 'rb')
genomic = pickle.load(f)
f.close()
y = pd.read_csv(labfn, header=0, index_col=0)

print('matrix shape: {0} x {1}'.format(*genomic.shape))

# read in demographic information as confounding variables
pts = pd.read_csv(ptsfn, index_col=0)

# make sure that the genetic data and confounding variables match each other
pts_sel = pts.reindex(index = genomic.index)
pts_sel.fillna(0, inplace=True)
sel_pts = np.array(pts_sel)

# read in train val test split, use pre-generated indices for reproducibility
train_indices = pd.read_csv('%s/train_indices_0.1val_0.2te.csv' % (indexdir), header=None) #  _5run
test_indices = pd.read_csv('%s/test_indices_0.1val_0.2te.csv' % (indexdir), header=None)
val_indices = pd.read_csv('%s/val_indices_0.1val_0.2te.csv' % (indexdir), header=None)
X = np.array(genomic)
y, yuniques = pd.factorize(y.values.ravel(), sort=True)



r = 0
if genfn=='/projects/b1122/saya/scanmap_data/gene_thres_thres010_abs.pik':
    ncs = range(100,501,100) # range(50,501,50)
elif genfn=='/projects/b1122/saya/scanmap_data/reg_thres.pik':
    ncs = range(60,301,60) # range(30,301,30)

X = X.astype(float)
device = torch.device(devstr)

train_index = train_indices[r]; val_index = val_indices[r]; test_index = test_indices[r]
X_train, X_val, X_test = X[train_index], X[val_index], X[test_index]
y_train, y_val, y_test = y[train_index], y[val_index], y[test_index]
pts_tr, pts_val, pts_te = sel_pts[train_index], sel_pts[val_index], sel_pts[test_index]    

print('nc,wcls,C,lr,tr bal acc,val bal acc,te acc,te bal acc,te precis,te recall,te f1,w2,b2,celoss') # ,best iter,mse,mse tr,mse val,mse te
for lr in [0.001, 0.01]: # 0.00001, 0.0001, 
    for nc in ncs:
        for C in [0.01, 0.1, 1, 10, 100]: # , 0.1, 1, 10, 100, 1000, 0.001, 
            for wcls in [0.01, 0.1, 1]: # , 2, 10, 50, 0.5, 
                #### to implement k-fold cross-validation, perform the splits here #### 
                #### run a for-loop here and record performance for every fold to average it later #### 
                fn = '%s/scanmap_k%d_wcls%s_C%s_lr%s.p' % (outdir, nc, wcls, C, lr) # scanmap%d/s%d/ # niter, seed, 
                m = ScanMap(
                    np.vstack((X_train, X_val, X_test)), 
                    cf = pts_tr, cfval = pts_val, y = y_train, yval = y_val, 
                    k=nc, n_iter=niter, weight_decay=0, lr=lr, wcls=wcls, C=C, seed=2*722019+seed, fn=fn, device=device
                ) # 2*722019+seed is just to have a large odd number for seeding that is recommended for generating random numbers, fixed for reproducibility
                [X_tr_nmf, X_val_nmf, X_te_nmf, H] = m.fit_transform()

                chkpt = torch.load(fn)
                m.load_state_dict(chkpt['state_dict'])
                # best_iter = chkpt['epoch']
                accval = chkpt['best_val_acc']
                m.eval()

                y_tr_pred = m.predict(X_tr_nmf, pts_tr)
                y_te_pred = m.predict(X_te_nmf, pts_te)
                acctr = balanced_accuracy_score(y_train, y_tr_pred)
                accte = accuracy_score(y_test, y_te_pred)
                balaccte = balanced_accuracy_score(y_test, y_te_pred)
                preciste = precision_score(y_test, y_te_pred)
                recallte = recall_score(y_test, y_te_pred)
                f1te = f1_score(y_test, y_te_pred)
                
                w2 = np.square(m.state_dict()['fc.weight'].cpu().numpy()).sum(axis=None)
                b2 = np.square(m.state_dict()['fc.bias'].cpu().numpy()).sum(axis=None)
                w = 1 / pd.Series(y_train).value_counts(normalize=True).sort_index().to_numpy()
                vce = chkpt['celoss']

                # err = np.vstack((X_train, X_val, X_test)) - np.vstack((X_tr_nmf, X_val_nmf, X_te_nmf)) @ H
                # err_tr = X_train - X_tr_nmf @ H
                # err_val = X_val - X_val_nmf @ H
                # err_te = X_test - X_te_nmf @ H
                # mse = np.square(err).mean(axis=None)
                # mse_tr = np.square(err_tr).mean(axis=None)
                # mse_val = np.square(err_val).mean(axis=None)
                # mse_te = np.square(err_te).mean(axis=None)
                
                print('%d,%s,%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f' % 
                    (nc, wcls, C, lr, acctr, accval, accte, balaccte, preciste, recallte, f1te, w2, b2, vce)) # best_iter, mse, mse_tr, mse_val, mse_te

# print("Ending run ... {}".format(datetime.datetime.now()))
