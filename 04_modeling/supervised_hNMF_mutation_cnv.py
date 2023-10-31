import numpy as np
import pandas as pd
import torch
from torch import nn
from torch import optim
from sklearn import metrics
import matplotlib.pyplot as plt

from suphNMF import suphNMF

X_mut = pd.read_csv('/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix/20230423_signature_results/sbs_96_original_per_sample.csv', index_col=0)
X_cnv = pd.read_csv('/projects/b1131/saya/bbcar/data/02b_cnv/inhouse_signatures/inhouse_sig_batcheffect_rm_combat.csv', index_col=0)
y = pd.read_csv('/projects/b1131/saya/bbcar/data/clinical/bbcar_label_studyid_from_gatk_filenames.csv', index_col=0)

indexdir = '/projects/b1131/saya/bbcar/train_test_splits'

train_ix = pd.read_csv(f'{indexdir}/train_index.txt', header=None).to_numpy().ravel()
test_ix = pd.read_csv(f'{indexdir}/test_index.txt', header=None).to_numpy().ravel()

X_mut_train, X_mut_test = X_mut.loc[train_ix,:], X_mut.loc[test_ix,:]
X_cnv_train, X_cnv_test = X_cnv.loc[train_ix,:], X_cnv.loc[test_ix,:]
y_train, y_test = y.loc[train_ix,:], y.loc[test_ix,:]

loss_type = ['l2', 'kl'][0]

test = suphNMF(X_mut_train, X_cnv_train, y_train, n_iter=10000, lr=1e-3, clf_weight=1.)
test.fit()

fig, axs = plt.subplots(2, 3, figsize=(12, 8))

axs[0,0].plot(test.loss_record)
axs[0,0].set_ylabel('Loss')
axs[0,0].set_xlabel('Epochs')
axs[0,0].set_title('Overall Loss')

axs[0,1].plot(test.reconloss1_record)
axs[0,1].set_xlabel('Epochs')
axs[0,1].set_ylabel('Loss')
axs[0,1].set_title('Reconstruction Loss (Mutation)')

axs[0,2].plot(test.reconloss2_record)
axs[0,2].set_xlabel('Epochs')
axs[0,2].set_ylabel('Loss')
axs[0,2].set_title('Reconstruction Loss (CNV)')

axs[1,0].plot(test.ortholoss_record)
axs[1,0].set_xlabel('Epochs')
axs[1,0].set_ylabel('Loss')
axs[1,0].set_title('Orthogonality loss')

axs[1,1].plot(test.clf_roc_record)
axs[1,1].set_xlabel('Epochs')
axs[1,1].set_ylabel('Classification performance')
axs[1,1].set_title('Classification ROC-AUC')

axs[1,2].plot(test.clf_loss_record)
axs[1,2].set_xlabel('Epochs')
axs[1,2].set_ylabel('Loss')
axs[1,2].set_title('Classification Loss')

plt.tight_layout()

fig.savefig('/projects/b1131/saya/bbcar/plots/test_suphNMF_learning_curves.png')
plt.close()
