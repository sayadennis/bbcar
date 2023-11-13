import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

dout = '/projects/b1131/saya/bbcar/plots/mutation'

###################
#### Load data ####
###################

## Load mutational features
sbs = pd.read_csv('/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix/20230423_signature_results/SBS96/Samples.txt', sep='\t', index_col=0).T
sbs.index = [int(x.split('_')[0]) for x in sbs.index]

dbs = pd.read_csv('/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix/20230423_signature_results/DBS78/Samples.txt', sep='\t', index_col=0).T
dbs.index = [int(x.split('_')[0]) for x in dbs.index]

indel = pd.read_csv('/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix/20230423_signature_results/ID83/Samples.txt', sep='\t', index_col=0).T
indel.index = [int(x.split('_')[0]) for x in indel.index]

data = { 
    'SBS' : sbs, 'DBS' : dbs, 'INDEL' : indel, 'ALL' : pd.concat((sbs, dbs, indel), axis=1)
}

## Load labels
label = pd.read_csv('/projects/b1131/saya/bbcar/data/clinical/bbcar_label_studyid_from_gatk_filenames.csv', index_col=0)

## Import sequencing source information
with open('/projects/b1131/saya/bbcar/data/sample_ids_uchicago.txt', 'r') as f:
    uchicago_samples = [int(x.strip()) for x in f.readlines()]

seqsource = pd.DataFrame(np.array([1 if x in uchicago_samples else 0 for x in sbs.index]).reshape(-1,1), index=sbs.index, columns=['label'])

## Load clinical features
clin = pd.read_csv('/projects/b1131/saya/bbcar/data/clinical/BBCaRDatabaseNU09B2_DATA_LABELS_2023-11-12_1525_utf8.csv', sep=',', index_col=None)
clin['study_id'] = [int(x.split('-')[0][3:]) for x in clin['Study ID']]
clin = clin.iloc[[clin.loc[i,'study_id'] in sbs.index for i in clin.index],:]

##################################################
#### Check data consistency and distributions ####
##################################################

# For a single study participant, is there a single label for ER, PR and HER2?
multiple_labels_er = []
multiple_labels_pr = []
multiple_labels_her2 = []
for i in sbs.index:
    subclin = clin.iloc[[x==i for x in clin.study_id.values],:]
    if len(subclin['ER Range'].dropna().unique())>1:
        multiple_labels_er.append(i)
    if len(subclin['PR Range'].dropna().unique())>1:
        multiple_labels_pr.append(i)
    if len(subclin['Her2 Range'].dropna().unique())>1:
        multiple_labels_her2.append(i)

if np.any([len(x)>0 for x in [multiple_labels_er, multiple_labels_pr, multiple_labels_her2]]):
    print('WARNING: Discrepancy in some of the biomarker labels from RedCAP.')

labels_all = pd.DataFrame(
    index=sbs.index,
    columns=['Case/Control', 'ER', 'PR', 'HER2', 'Sequenced at']
)

for i in sbs.index:
    labels_all.loc[i,'Case/Control'] = label.loc[i,'label'] if i in label.index else None
    labels_all.loc[i,'ER'] = clin.iloc[clin.study_id.values==i,:]['ER Range'].dropna().values[0] if len(clin.iloc[clin.study_id.values==i,:]['ER Range'].dropna().values)>0 else None
    labels_all.loc[i,'PR'] = clin.iloc[clin.study_id.values==i,:]['PR Range'].dropna().values[0] if len(clin.iloc[clin.study_id.values==i,:]['PR Range'].dropna().values)>0 else None
    labels_all.loc[i,'HER2'] = clin.iloc[clin.study_id.values==i,:]['Her2 Range'].dropna().values[0] if len(clin.iloc[clin.study_id.values==i,:]['Her2 Range'].dropna().values)>0 else None

labels_all['Sequenced at'] = seqsource.values.ravel()

# Make sure that no controls have biomarker labels and check distributions across categories
labels_all.groupby(['Case/Control', 'ER', 'PR', 'HER2']).size()

#####################################################
#### Plot PCA to ensure there is no batch effect ####
#####################################################

fig, ax = plt.subplots(1, 4, figsize=(12,3))

for i, key in enumerate(data.keys()):
    # Run PCA
    pca = PCA(n_components=2)
    X_tf = pca.fit_transform(data[key])
    # Plot
    ax[i].scatter(X_tf[:,0][labels_all['Sequenced at'].values==0], X_tf[:,1][labels_all['Sequenced at'].values==0], label='Indiana', alpha=0.5)
    ax[i].scatter(X_tf[:,0][labels_all['Sequenced at'].values==1], X_tf[:,1][labels_all['Sequenced at'].values==1], label='UChicago', alpha=0.5)
    ax[i].set_ylabel('PC2')
    ax[i].set_title(f'{key}')
    if key == 'SBS':
        ax[i].legend()
    if i==len(data.keys()):
        ax[i].xlabel('PC1')

fig.suptitle('PCA colored by sequencing source')
plt.tight_layout()
fig.savefig(f'{dout}/mutsig_features_pca_seq_source.png')
plt.close()

############################################
#### Plot PCA with corresponding labels ####
############################################

for i in labels_all.index:
    if labels_all.loc[i,'Case/Control']==0:
        labels_all.loc[i,'ER'] = 'Control'
        labels_all.loc[i,'PR'] = 'Control'
        labels_all.loc[i,'HER2'] = 'Control'

label_ix_map = {
    'ER' : {'Control': 0, '0% Negative': 1, '1-9% Low Positive': 2, '10% Positive': 3},
    'PR' : {'Control': 0, '0% Negative': 1, '1-9% Low Positive': 2, '10% Positive': 3},
    'HER2' : {'Control': 0, '0 Negative': 1, '1+ Negative': 2, '2+ Equivocal': 3, '3+ Positive' : 4},
}

colors = ['mediumseagreen', 'cornflowerblue', 'darkviolet', 'fuchsia', 'red']

fig, ax = plt.subplots(3, 4, figsize=(12,9))

for j, key in enumerate(data.keys()):
    # Run PCA
    pca = PCA(n_components=2)
    X_tf = pca.fit_transform(data[key])
    # Plot
    for i, coloring in enumerate(['ER', 'PR', 'HER2']):
        # Get string to integer mapping
        for k, v in label_ix_map[coloring].items():
            ax[i,j].scatter(
                X_tf[:,0][labels_all[coloring].map(label_ix_map[coloring]).values==v], 
                X_tf[:,1][labels_all[coloring].map(label_ix_map[coloring]).values==v], 
                label=k, alpha=0.5, s=10, c=colors[v])
            ax[i,j].set_xticklabels([])
            ax[i,j].set_yticklabels([])
        if i==0:
            ax[i,j].set_title(f'{key}', fontsize=14)
        if i==2:
            ax[i,j].set_xlabel('PC1')
        if j==0: 
            ax[i,j].legend()
            ax[i,j].set_ylabel('PC2')

fig.suptitle('PCA colored by biomarker levels', fontsize=16)

fig.subplots_adjust(left=0.1)
plt.text(0.03, 0.68, 'Estrogen Receptor', color='black', fontsize=14, rotation='vertical', transform=plt.gcf().transFigure)
plt.text(0.03, 0.38, 'Progesterone Receptor', color='black', fontsize=14, rotation='vertical', transform=plt.gcf().transFigure)
plt.text(0.03, 0.20, 'HER2', color='black', fontsize=14, rotation='vertical', transform=plt.gcf().transFigure)

fig.savefig(f'{dout}/mutsig_features_pca_biomarker.png')
plt.close()


