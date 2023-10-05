import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

din = '/projects/b1131/saya/bbcar/data/02b_cnv/inhouse_signatures'
dout = '/projects/b1131/saya/bbcar/plots/cnv'

################################################
#### Original data w/o batch effect removal ####
################################################

data = pd.read_csv(f'{din}/seglen_ampdel_category_call_counts_per_sample.csv', index_col=0)

pca = PCA(n_components=10)
data_tf = pca.fit_transform(data.values)

## Get the sample IDs of U Chicago samples 
with open('/projects/b1131/saya/bbcar/data/sample_ids_uchicago.txt', 'r') as f:
    uchicago_samples = [int(x.strip()) for x in f.readlines()]

#### Plot explained variance ####

fig, ax = plt.subplots()
bars = ax.bar(np.arange(10), pca.explained_variance_ratio_[:10])

for bar in bars:
    height = bar.get_height()
    ax.annotate(f'{height:.1e}', 
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),  # 3 points vertical offset
                textcoords="offset points", rotation=30,
                ha='left', va='bottom')

ax.set_ylabel('Explained variance ratio')
ax.set_xlabel('Principal component')
ax.set_ylim(0,0.97)
ax.set_xticks(np.arange(10))
ax.set_xticklabels(np.arange(1,11))
ax.set_title('Explained variance ratios of in-house CNV signatures')
plt.tight_layout()
fig.savefig(f'{dout}/cnv_inhouse_signature_pca_explained_variance_ratio.png')
plt.close()

#### Plot scatterplot of PCs ####

labels_num = [1 if x in uchicago_samples else 0 for x in data.index]

pc1 = data_tf[:,0]
pc2 = data_tf[:,1]

fig, ax = plt.subplots()
ax.scatter(pc1[np.array(labels_num)==0], pc2[np.array(labels_num)==0], label='Indiana')
ax.scatter(pc1[np.array(labels_num)==1], pc2[np.array(labels_num)==1], label='U Chicago')
ax.legend(title='Source') #  *scatter.legend_elements(), title='Source')
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_title('PCA of in-house CN sig before batch effect removal')
fig.savefig(f'{dout}/inhouse_cnv_signature_pc1_and_pc2_scatter.png')
plt.close()

pc3 = data_tf[:,2]
pc4 = data_tf[:,3]

fig, ax = plt.subplots()
scatter = ax.scatter(pc3, pc4, c=labels_num, label=['Indiana', 'UChicago'])
ax.legend(*scatter.legend_elements(), title='Source')
ax.set_xlabel('PC3')
ax.set_ylabel('PC4')
ax.set_title('PCA of in-house CN sig before batch effect removal')
plt.savefig(f'{dout}/inhouse_cnv_signature_pc3_and_pc4_scatter.png')
plt.close()

############################################
#### Data with PCA batch effect removal ####
############################################

data = pd.read_csv(f'{din}/inhouse_sig_batcheffect_rm_pc1.csv', index_col=0)

pca = PCA(n_components=10)
data_tf = pca.fit_transform(data.values)

## Get the sample IDs of U Chicago samples 
with open('/projects/b1131/saya/bbcar/data/sample_ids_uchicago.txt', 'r') as f:
    uchicago_samples = [int(x.strip()) for x in f.readlines()]

#### Plot explained variance ####

fig, ax = plt.subplots()
bars = ax.bar(np.arange(10), pca.explained_variance_ratio_[:10])

for bar in bars:
    height = bar.get_height()
    ax.annotate(f'{height:.1e}', 
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),  # 3 points vertical offset
                textcoords="offset points", rotation=30,
                ha='left', va='bottom')

ax.set_ylabel('Explained variance ratio')
ax.set_xlabel('Principal component')
ax.set_ylim(0,0.97)
ax.set_xticks(np.arange(10))
ax.set_xticklabels(np.arange(1,11))
ax.set_title('Explained variance ratios of in-house CNV signatures post-PCA batch effect removal')
plt.tight_layout()
fig.savefig(f'{dout}/cnv_inhouse_signature_pca_explained_variance_ratio_batchrm_pc1.png')
plt.close()

#### Plot scatterplot of PCs ####

labels_num = [1 if x in uchicago_samples else 0 for x in data.index]

pc1 = data_tf[:,0]
pc2 = data_tf[:,1]

fig, ax = plt.subplots()
ax.scatter(pc1[np.array(labels_num)==0], pc2[np.array(labels_num)==0], label='Indiana')
ax.scatter(pc1[np.array(labels_num)==1], pc2[np.array(labels_num)==1], label='U Chicago')
ax.legend(title='Source')
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_title('PCA of in-house CN sig post-PCA batch effect removal')
fig.savefig(f'{dout}/inhouse_cnv_signature_pc1_and_pc2_scatter_batchrm_pc1.png')
plt.close()

###############################################
#### Data with ComBat batch effect removal ####
###############################################

data = pd.read_csv(f'{din}/inhouse_sig_batcheffect_rm_combat.csv', index_col=0)

pca = PCA(n_components=10)
data_tf = pca.fit_transform(data.values)

## Get the sample IDs of U Chicago samples 
with open('/projects/b1131/saya/bbcar/data/sample_ids_uchicago.txt', 'r') as f:
    uchicago_samples = [int(x.strip()) for x in f.readlines()]

#### Plot explained variance ####

fig, ax = plt.subplots()
bars = ax.bar(np.arange(10), pca.explained_variance_ratio_[:10])

for bar in bars:
    height = bar.get_height()
    ax.annotate(f'{height:.1e}', 
                xy=(bar.get_x() + bar.get_width() / 2, height),
                xytext=(0, 3),  # 3 points vertical offset
                textcoords="offset points", rotation=30,
                ha='left', va='bottom')

ax.set_ylabel('Explained variance ratio')
ax.set_xlabel('Principal component')
ax.set_ylim(0,0.97)
ax.set_xticks(np.arange(10))
ax.set_xticklabels(np.arange(1,11))
ax.set_title('Explained variance ratios of in-house CNV signatures post-combat batch effect removal')
plt.tight_layout()
fig.savefig(f'{dout}/cnv_inhouse_signature_pca_explained_variance_ratio_batchrm_combat.png')
plt.close()

#### Plot scatterplot of PCs ####

labels_num = [1 if x in uchicago_samples else 0 for x in data.index]

pc1 = data_tf[:,0]
pc2 = data_tf[:,1]

fig, ax = plt.subplots()
ax.scatter(pc1[np.array(labels_num)==0], pc2[np.array(labels_num)==0], label='Indiana')
ax.scatter(pc1[np.array(labels_num)==1], pc2[np.array(labels_num)==1], label='U Chicago')
ax.legend(title='Source')
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_title('PCA of in-house CN sig post-combat batch effect removal')
fig.savefig(f'{dout}/inhouse_cnv_signature_pc1_and_pc2_scatter_batchrm_combat.png')
plt.close()


