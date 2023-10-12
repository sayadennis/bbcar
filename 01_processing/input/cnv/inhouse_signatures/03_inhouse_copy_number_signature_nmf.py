import numpy as np
import pandas as pd
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import silhouette_score
from sklearn.decomposition import NMF
import matplotlib.pyplot as plt

#############################
#### Load feature matrix ####
#############################

dn = '/projects/b1131/saya/bbcar/data/02b_cnv/inhouse_signatures'
plot_dn = '/projects/b1131/saya/bbcar/plots/cnv'
fn = 'inhouse_cn_features_batcheffect_rm_combat.csv'

X = pd.read_csv(f'{dn}/{fn}', index_col=0)

# EXPERIMENTAL - make all non-negative (since matrix contains negative values following batch effect removal)
X -= np.min(X.values)

#############################################################
#### Perform NMF and determine best number of components ####
#############################################################

k_list = np.arange(2,16)
nmf_scores = pd.DataFrame(index=k_list, columns=['recon_error', 'stability'])

for k in k_list:
    nmf = NMF(n_components=k, init='random', random_state=0, max_iter=2000)
    W = nmf.fit_transform(X)
    H = nmf.components_
    X_recon = np.dot(W, H)
    recon_error = np.linalg.norm(np.array(X - X_recon), ord='fro') # frobenious norm as described in the COSMIC paper
    stability = silhouette_score(X, np.argmax(W, axis=1)) # stability is measured by silhouette score - cluster assignment determined by W matrix
    nmf_scores.loc[k,'recon_error'] = recon_error
    nmf_scores.loc[k,'stability'] = stability

scaler = MinMaxScaler()
scaler.fit(nmf_scores)
scaled_scores = pd.DataFrame(scaler.transform(nmf_scores), index=nmf_scores.index, columns=nmf_scores.columns)

best_k = (scaled_scores['stability'] - scaled_scores['recon_error']).idxmax() # index that gives the max difference

fig, ax1 = plt.subplots()

color = 'tab:red'
ax1.set_xlabel('Number of components')
ax1.set_ylabel('Stability', color=color)
ax1.plot(k_list, nmf_scores['stability'], color=color, marker='s')
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx() # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('Reconstruction error', color=color)
ax2.plot(k_list, nmf_scores['recon_error'], color=color, marker='s')
ax2.tick_params(axis='y', labelcolor=color)

ax2.axvline(best_k, linestyle='--', color='orange')

fig.suptitle('Stability and Reconstruction Error of NMF')

fig.tight_layout()  # otherwise the right y-label is slightly clipped
fig.savefig(f'{plot_dn}/inhouse_cn_signature_stability_recon_plot.png')
plt.close()

#############################
#### Save feature matrix ####
#############################

nmf = NMF(n_components=best_k, init='random', random_state=0, max_iter=2000)
W = nmf.fit_transform(X)
H = nmf.components_

pd.DataFrame(W, index=X.index, columns=np.arange(W.shape[1])).to_csv(f'{dn}/inhouse_cnv_sig_per_sample.csv', header=True, index=True)
pd.DataFrame(H, index=np.arange(H.shape[0]), columns=X.columns).T.to_csv(f'{dn}/inhouse_cnv_signature_features.csv', header=True, index=True)
