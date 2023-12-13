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
cv_results, best_params = nmf.crossval_fit(
    n_iters=[1000,2000], lrs=[1e-2], clf_weights=[1e-1,1e+0,1e+1,1e+2],
    #weight_decays=[1e-3], 
    #ortho_weights=[1.],
)

eval_metrics = nmf.evaluate(
    torch.tensor(X_mut_test.values, dtype=torch.float32),
    torch.tensor(X_cnv_test.values, dtype=torch.float32),
    torch.tensor(y_test.values, dtype=torch.float32),
)

y_pred = eval_metrics['y_pred'].detach().numpy()
y_score = eval_metrics['y_score'].detach().numpy()

print(f'Training ROC-AUC: {nmf.train_roc}')
print(f'Cross-validation ROC-AUC: {nmf.cv_mean_roc}')
print(f'Test ROC-AUC: {metrics.roc_auc_score(y_test, y_score)}')
print(f'Test precision: {metrics.precision_score(y_test, y_pred)}')
print(f'Test recall: {metrics.recall_score(y_test, y_pred)}')
print(f'Test F1: {metrics.f1_score(y_test, y_pred)}')

with open('/projects/b1131/saya/bbcar/model_interpretations/test_cv_results.p', 'wb') as f:
    pickle.dump(cv_results, f)

with open('/projects/b1131/saya/bbcar/model_interpretations/test_best_params.p', 'wb') as f:
    pickle.dump(best_params, f)

with open('/projects/b1131/saya/bbcar/model_interpretations/test_hNMF.p', 'wb') as f:
    pickle.dump(nmf, f)

with open('/projects/b1131/saya/bbcar/model_interpretations/test_eval_metrics.p', 'wb') as f:
    pickle.dump(eval_metrics, f)

