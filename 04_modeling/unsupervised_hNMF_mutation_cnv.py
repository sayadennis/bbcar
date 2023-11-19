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

#loss_type = ['l2', 'kl'][i]

def plot_learning(weight_decay, clf_weight, ortho_weight):
    test = suphNMF(X_mut_train, X_cnv_train, n_iter=10000, lr=1e-3, weight_decay=weight_decay, clf_weight=clf_weight, ortho_weight=ortho_weight)
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
    
    for i in range(2):
        for j in range(3):
            axs[i,j].set_xlabel('Epochs')
            if not (i==1) & (j==1):
                axs[i,j].set_ylabel('Loss')

    fig.suptitle(f'Learning curves - weight decay {weight_decay}; classification weight {clf_weight}; orthogonality weight {ortho_weight}')
    plt.tight_layout()
    fig.savefig(f'{plotdir}/test_hNMF_learning_curves_decay{weight_decay}_clf{clf_weight}_ortho{ortho_weight}.png')
    plt.close()

clf_weight = 1.

for weight_decay in [1e-4, 1e-3, 1e-2]:
    for ortho_weight in [1e-1, 1e+0, 1e+1, 1e+2]:
        plot_learning(weight_decay, clf_weight, ortho_weight)

k_list = np.arange(2,25)
nmf_scores = pd.DataFrame(index=k_list, columns=['recon_error1', 'recon_error2', 'stability1', 'stability2'])
for k in k_list:
    nmf = suphNMF(
        X_mut_train, X_cnv_train, # y_train, 
        n_components=k, n_iter=int(1e+4), weight_decay=1e-3, 
        clf_weight=1e+0, ortho_weight=1e+0
    )
    nmf.fit()
    W_train = nmf.W.detach()
    X1_train_recon = torch.mm(W_train, nmf.H1.detach())
    X2_train_recon = torch.mm(W_train, nmf.H2.detach())
    # Calculate reconstruction error 
    #recon_error1 = np.linalg.norm(np.array(nmf.X1 - X1_train_recon), ord='fro') # Frobenious norm or Cosine?!
    #recon_error2 = np.linalg.norm(np.array(nmf.X2 - X2_train_recon), ord='fro') # Frobenious norm or Cosine?!
    recon_error1 = np.mean(np.diag(cdist(X_mut_train, X1_train_recon, metric='cosine')))
    recon_error2 = np.mean(np.diag(cdist(X_cnv_train, X2_train_recon, metric='cosine')))
    # Calculate stability 
    stability1 = metrics.silhouette_score(nmf.X1, np.argmax(W_train, axis=1)) # stability is measured by silhouette score - cluster assignment determined by W matrix
    stability2 = metrics.silhouette_score(nmf.X2, np.argmax(W_train, axis=1))
    # Record values
    nmf_scores.loc[k,'recon_error1'] = recon_error1
    nmf_scores.loc[k,'recon_error2'] = recon_error2
    nmf_scores.loc[k,'stability1'] = stability1
    nmf_scores.loc[k,'stability2'] = stability2

###############################
#### Determine the best K  ####
###############################

scaler = MinMaxScaler()
scaler.fit(nmf_scores)
scaled_scores = pd.DataFrame(scaler.transform(nmf_scores), index=nmf_scores.index, columns=nmf_scores.columns)

best_k_1 = (scaled_scores['stability1'] - scaled_scores['recon_error1']).idxmax() # index that gives the max difference
best_k_2 = (scaled_scores['stability2'] - scaled_scores['recon_error2']).idxmax() # index that gives the max difference

ix_to_name = {1: 'Mutation', 2: 'CNV'}
ix_to_bestk = {1: best_k_1, 2: best_k_2}

fig, axs1 = plt.subplots(1, 2, figsize=(15,5))

for i in range(2):
    color = 'tab:red'
    axs1[i].set_xlabel('Number of components')
    axs1[i].set_ylabel('Stability', color=color)
    axs1[i].plot(k_list, nmf_scores[f'stability{i+1}'], color=color, marker='s')
    axs1[i].tick_params(axis='y', labelcolor=color)
    # instantiate second axes that shares the same x-axis
    ax2 = axs1[i].twinx()
    # plot on the second axes
    color = 'tab:blue'
    ax2.set_ylabel('Reconstruction error', color=color)
    ax2.plot(k_list, nmf_scores[f'recon_error{i+1}'], color=color, marker='s')
    ax2.tick_params(axis='y', labelcolor=color)
    # add vertical line on the best K
    ax2.axvline(ix_to_bestk[i+1], linestyle='--', color='orange')
    # set title
    axs1[i].set_title(f'{ix_to_name[i+1]}')

fig.suptitle('Stability and Reconstruction Error of hybrid NMF')

fig.tight_layout()  # otherwise the right y-label is slightly clipped
fig.savefig(f'{plotdir}/hNMF_stability_recon_plot.png')
plt.close()

###############################
#### Save best model and W #### # TODO: Needs re-writing!! 
###############################

best_k = 4 # visually determined from mutation and CNV plot (better ways?)

nmf = suphNMF(
    X_mut_train, X_cnv_train, # y_train, 
    n_components=best_k, n_iter=int(1e+4), weight_decay=1e-3, 
    clf_weight=1e+0, ortho_weight=1e+0
)   
nmf.fit()
W_train = nmf.W.detach()

W_test = nmf.transform(X_mut_test, X_cnv_test, y_test).detach().numpy()
W_train = nmf.W.detach().numpy()

W = pd.concat(
    (pd.DataFrame(W_train, index=X_mut_train.index),
    pd.DataFrame(W_test, index=X_mut_test.index)),
    axis=0
)

W.to_csv(f'{datadir}/combined_mutation_cnv/hNMF_learned_W.csv', index=True, header=True)

with open(f'/projects/b1131/saya/bbcar/model_interpretations/suphNMF/hNMF.p', 'wb') as f:
    pickle.dump(nmf, f)

# Calculate ROC-AUC
#y_test_pred = torch.sigmoid(nmf.predict(nmf.transform(X_mut_test, X_cnv_test, y_test)))
#print(f'ROC-AUC: ', metrics.roc_auc_score(y_test, y_test_pred.detach().numpy()))
