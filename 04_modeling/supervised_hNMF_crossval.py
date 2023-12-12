import pickle

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch import optim
from scipy.spatial.distance import cdist
from sklearn import metrics
from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt

from suphNMF import suphNMF

datadir = '/projects/b1131/saya/bbcar/data'
plotdir = '/projects/b1131/saya/bbcar/plots/suphNMF_learning_curves'

X_mut = pd.read_csv(
    f'{datadir}/02a_mutation/08_feature_matrix/20230423_signature_results/sbs_96_original_per_sample.csv', 
    index_col=0
)
X_cnv = pd.read_csv(
    f'{datadir}/02b_cnv/inhouse_signatures/inhouse_sig_batcheffect_rm_combat.csv', 
    index_col=0
)
y = pd.read_csv(
    f'{datadir}/clinical/bbcar_label_studyid_from_gatk_filenames.csv', 
    index_col=0
)

indexdir = '/projects/b1131/saya/bbcar/train_test_splits'

train_ix = pd.read_csv(f'{indexdir}/train_index.txt', header=None).to_numpy().ravel()
test_ix = pd.read_csv(f'{indexdir}/test_index.txt', header=None).to_numpy().ravel()

X_mut_train, X_mut_test = X_mut.loc[train_ix,:], X_mut.loc[test_ix,:]
X_cnv_train, X_cnv_test = X_cnv.loc[train_ix,:], X_cnv.loc[test_ix,:]
y_train, y_test = y.loc[train_ix,:], y.loc[test_ix,:]

nmf = suphNMF(X_mut_train, X_cnv_train, y_train, n_components=7)
cv_results, best_params = nmf.fit_cv()

with open('/projects/b1131/saya/bbcar/model_interpretations/test_cv_results.p', 'wb') as f:
    pickle.dump(cv_results, f)

with open('/projects/b1131/saya/bbcar/model_interpretations/test_best_params.p', 'wb') as f:
    pickle.dump(best_params, f)

with open('/projects/b1131/saya/bbcar/model_interpretations/test_hNMF.p', 'wb') as f:
    pickle.dump(nmf, f)

