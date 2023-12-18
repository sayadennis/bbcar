import os
import pickle

import numpy as np
import pandas as pd
import torch
from torch import nn
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
) * 10
X_cnv = pd.read_csv(
    f'{datadir}/02b_cnv/inhouse_signatures/inhouse_sig_batcheffect_rm_combat.csv', 
    index_col=0
) * 10
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

dout = '/projects/b1131/saya/bbcar/model_interpretations/supervised_hNMF'

with open(f'{dout}/suphNMF.p', 'wb') as f:
    pickle.dump(nmf, f)

with open(f'{dout}/eval_metrics.p', 'wb') as f:
    pickle.dump(eval_metrics, f)

## Get confusion matrix
cf_mx = pd.DataFrame(
    metrics.confusion_matrix(y_test, eval_metrics['y_pred'].detach().numpy()),
    index = ['True control', 'True case'],
    columns = ['Predicted control', 'Predicted case']
)

## Plot ROC
fpr, tpr, _ = metrics.roc_curve(y_test, eval_metrics['y_score'].detach().numpy())
auc = metrics.auc(fpr, tpr)

fig, ax = plt.subplots()
ax.plot(fpr, tpr, color="darkorange", lw=2, label="ROC curve (area = %0.3f)" % auc,)
ax.plot([0, 1], [0, 1], color="navy", lw=2, linestyle="--")
ax.set_xlim([0.0, 1.0])
ax.set_ylim([0.0, 1.05])
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
ax.set_title("Receiver operating characteristic on test set")
ax.legend(loc="lower right")

fig.savefig(f'{dout}/suphNMF_roc_curve_test_set.png')
plt.close()
