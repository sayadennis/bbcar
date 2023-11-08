import pickle

import numpy as np
import pandas as pd
import torch
from torch import nn
from torch import optim
from sklearn import metrics
import matplotlib.pyplot as plt

from suphNMF import suphNMF

datadir = '/projects/b1131/saya/bbcar/data'

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

loss_type = ['l2', 'kl'][0]

def plot_learning(weight_decay, clf_weight, ortho_weight):
    test = suphNMF(X_mut_train, X_cnv_train, y_train, n_iter=10000, lr=1e-3, weight_decay=weight_decay, clf_weight=clf_weight, ortho_weight=ortho_weight)
    test.fit()
    
    fig, axs = plt.subplots(2, 3, figsize=(12, 8))
    
    axs[0,0].plot(test.loss_record)
    axs[0,0].set_title('Overall Loss')
    
    axs[0,1].plot(test.reconloss1_record)
    axs[0,1].set_title('Reconstruction Loss (Mutation)')
    
    axs[0,2].plot(test.reconloss2_record)
    axs[0,2].set_title('Reconstruction Loss (CNV)')
    
    axs[1,0].plot(test.ortholoss_record)
    axs[1,0].set_title('Orthogonality loss')
    
    axs[1,1].plot(test.clf_roc_record)
    axs[1,1].set_ylabel('Classification performance')
    axs[1,1].set_title('Classification ROC-AUC')
    
    axs[1,2].plot(test.clf_loss_record)
    axs[1,2].set_title('Classification Loss')
    
    for i in range(3):
        for j in range(2):
            axs[i,j].set_xlabel('Epochs')
            if not (i==1) & (j==1):
                axs[i,j].set_ylabel('Loss')

    fig.suptitle(f'Learning curves - weight decay {weight_decay}; classification weight {clf_weight}; orthogonality weight {ortho_weight}')
    plt.tight_layout()
    fig.savefig(f'/projects/b1131/saya/bbcar/plots/suphNMF_learning_curves/test_suphNMF_learning_curves_decay{weight_decay}_clf{clf_weight}_ortho{ortho_weight}.png')
    plt.close()

#for weight_decay in [1e-3, 1e-1, 1e+1]:
#    for clf_weight in [1e+0, 1e+2, 1e+3, 1e+4]:
#        for ortho_weight in [1e+0, 1e+2, 1e+3, 1e+4]:
#            plot_learning(weight_decay, clf_weight, ortho_weight)

test = suphNMF(X_mut_train, X_cnv_train, y_train, n_iter=int(1e+4), lr=1e-3, weight_decay=.1, clf_weight=1e+3, ortho_weight=1e+3)
test.fit()

W_test = test.transform(X_mut_test, X_cnv_test, y_test).detach().numpy()
W_train = test.W.detach().numpy()

W = pd.concat(
    (pd.DataFrame(W_train, index=X_mut_train.index),
    pd.DataFrame(W_test, index=X_mut_test.index)),
    axis=0
)

W.to_csv(f'{datadir}/combined_mutation_cnv/test_learned_W.csv', index=True, header=True)

with open(f'/projects/b1131/saya/bbcar/model_interpretations/suphNMF/suphNMF.p', 'wb') as f:
    pickle.dump(test, f)

# Calculate ROC-AUC
#y_test_pred = torch.sigmoid(test.predict(test.transform(X_mut_test, X_cnv_test, y_test)))
#print(f'ROC-AUC: ', metrics.roc_auc_score(y_test, y_test_pred.detach().numpy()))
