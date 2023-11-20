import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
import seaborn as sns

###################
#### Load data ####
###################

din = '/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix/20230423_signature_results'

## Load genomic features
data = {
    'De Novo' : {},
    'COSMIC' : {},
}

artifact_sigs = [f'SBS{i}' for i in [27,43,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,95]] + [f'DBS{i}' for i in [14]]

for sigcat in ['SBS96', 'DBS78', 'ID83']:
    samples = pd.read_csv(f'{din}/{sigcat}/Samples.txt', sep='\t', index_col=0).T
    denovo = pd.read_csv(f'{din}/{sigcat}/Suggested_Solution/{sigcat}_De-Novo_Solution/Signatures/{sigcat}_De-Novo_Signatures.txt', sep='\t', index_col=0)
    denovo_exposure = pd.DataFrame(np.dot(samples.values, denovo.values), index=[int(x.split('_')[0]) for x in samples.index], columns=denovo.columns)
    data['De Novo'][sigcat] = denovo_exposure
    data['COSMIC'][sigcat] = pd.read_csv(f'/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix/cosmic_signature_per_sample_{sigcat}.csv', index_col=0)
    data['COSMIC'][sigcat] = data['COSMIC'][sigcat].drop(artifact_sigs, axis=1, errors='ignore')

## Load clinical features
label = pd.read_csv('/projects/b1131/saya/bbcar/data/clinical/bbcar_label_studyid_from_gatk_filenames.csv', index_col=0)

clin = pd.read_csv('/projects/b1131/saya/bbcar/data/clinical/BBCaRDatabaseNU09B2_DATA_LABELS_2023-11-12_1525_utf8.csv', sep=',', index_col=None)
clin['study_id'] = [int(x.split('-')[0][3:]) for x in clin['Study ID']]
clin = clin.iloc[[clin.loc[i,'study_id'] in data['De Novo']['SBS96'].index for i in clin.index],:]

labels_all = pd.DataFrame(
    index=data['De Novo']['SBS96'].index,
    columns=['Case/Control', 'ER', 'PR', 'HER2', 'Sequenced at']
)

for i in data['De Novo']['SBS96'].index:
    labels_all.loc[i,'Case/Control'] = label.loc[i,'label'] if i in label.index else None
    labels_all.loc[i,'ER'] = clin.iloc[clin.study_id.values==i,:]['ER Range'].dropna().values[0] if len(clin.iloc[clin.study_id.values==i,:]['ER Range'].dropna().values)>0 else None
    labels_all.loc[i,'PR'] = clin.iloc[clin.study_id.values==i,:]['PR Range'].dropna().values[0] if len(clin.iloc[clin.study_id.values==i,:]['PR Range'].dropna().values)>0 else None
    labels_all.loc[i,'HER2'] = clin.iloc[clin.study_id.values==i,:]['Her2 Range'].dropna().values[0] if len(clin.iloc[clin.study_id.values==i,:]['Her2 Range'].dropna().values)>0 else None

for i in labels_all.index:
    if labels_all.loc[i,'Case/Control']==0:
        labels_all.loc[i,'ER'] = 'Control'
        labels_all.loc[i,'PR'] = 'Control'
        labels_all.loc[i,'HER2'] = 'Control'

###########################
#### Plot violin plots ####
###########################

sig_type = 'COSMIC'

sbs = data[sig_type]['SBS96']

fig, ax = plt.subplots(figsize=(16,4))

pos = np.arange(sbs.shape[1] * 4)
pos = pos[pos % 4 != 0]  # define positions to plot (this will look like [1,2,3,5,6,7,...])

plot_data = []

for signame in sbs.columns:
    for bio_group in ['Control', '0% Negative', '10% Positive']:
        plot_data.append(sbs.iloc[labels_all['ER'].values==bio_group,:][signame].values)

violins = ax.violinplot(
    plot_data, pos, widths=0.7, showmedians=True, showextrema=True
)

for i in [x for x in range(len(pos)) if x % 3 == 0]:
    violins['bodies'][i].set_facecolor('blue')

for i in [x for x in range(len(pos)) if x % 3 == 1]:
    violins['bodies'][i].set_facecolor('green')

for i in [x for x in range(len(pos)) if x % 3 == 2]:
    violins['bodies'][i].set_facecolor('orange')

ax.set_xticks([4*x+2 for x in np.arange(sbs.shape[1])])
ax.set_xticklabels(sbs.columns, ha='right', rotation=30)
ax.set_ylabel('Signature Exposure Levels')
ax.set_title('SBS Mutational Signature Exposure by ER status')
ax.legend(['Controls', 'ER Positive', 'ER Negative'], loc='upper right')

plt.tight_layout()
fig.savefig(f'/projects/b1131/saya/bbcar/plots/mutation/violin_sbs_exposure_by_er_status.png')
plt.close()

# Quick statistical test to see if ER+ seems to be more similar to controls (preliminary)

results = {}
for i, colname in enumerate(sbs.columns):
    results[colname] = {}
    t, p = ttest_ind(plot_data[i*3], plot_data[i*3+1]) # controls vs. ER+
    results[colname]['Control vs. ER+'] = {'t': t, 'p' : p}
    t, p = ttest_ind(plot_data[i*3+1], plot_data[i*3+2]) # ER+ vs. ER-
    results[colname]['ER+ vs. ER-'] = {'t' : t, 'p' : p}
    t, p = ttest_ind(plot_data[i*3], plot_data[i*3+2]) # Controls vs. ER-
    results[colname]['Control vs. ER-'] = {'t' : t, 'p' : p}

for signame in results.keys():
    print(f'#### {signame} ####')
    for comparison in results[signame].keys():
            print([f'{stat}={results[signame][comparison][stat]:.2f}' for stat in ['t', 'p']])

## DBS

dbs = data[sig_type]['DBS78']

fig, ax = plt.subplots(figsize=(3,4))

pos = np.arange(dbs.shape[1] * 4)
pos = pos[pos % 4 != 0]  # define positions to plot (this will look like [1,2,3,5,6,7,...])

plot_data = []

for signame in dbs.columns:
    for bio_group in ['Control', '0% Negative', '10% Positive']:
        plot_data.append(dbs.iloc[labels_all['ER'].values==bio_group,:][signame].values)

violins = ax.violinplot(
    plot_data, pos, widths=0.7, showmedians=True, showextrema=True
)

for i in [x for x in range(len(pos)) if x % 3 == 0]:
    violins['bodies'][i].set_facecolor('blue')

for i in [x for x in range(len(pos)) if x % 3 == 1]:
    violins['bodies'][i].set_facecolor('green')

for i in [x for x in range(len(pos)) if x % 3 == 2]:
    violins['bodies'][i].set_facecolor('orange')

ax.set_xticks([4*x+2 for x in np.arange(dbs.shape[1])])
ax.set_xticklabels(dbs.columns, ha='right', rotation=30)
ax.set_ylabel('Signature Exposure Levels')
ax.set_title('DBS Exposure')
#ax.legend(['Controls', 'ER Positive', 'ER Negative'], loc='upper right')

plt.tight_layout()
fig.savefig(f'/projects/b1131/saya/bbcar/plots/mutation/violin_dbs_exposure_by_er_status.png')
plt.close()

## ID83

indel = data[sig_type]['ID83']

fig, ax = plt.subplots(figsize=(10,4))

pos = np.arange(indel.shape[1] * 4)
pos = pos[pos % 4 != 0]  # define positions to plot (this will look like [1,2,3,5,6,7,...])

plot_data = []

for signame in indel.columns:
    for bio_group in ['Control', '0% Negative', '10% Positive']:
        plot_data.append(indel.iloc[labels_all['ER'].values==bio_group,:][signame].values)

violins = ax.violinplot(
    plot_data, pos, widths=0.7, showmedians=True, showextrema=True
)

for i in [x for x in range(len(pos)) if x % 3 == 0]:
    violins['bodies'][i].set_facecolor('blue')

for i in [x for x in range(len(pos)) if x % 3 == 1]:
    violins['bodies'][i].set_facecolor('green')

for i in [x for x in range(len(pos)) if x % 3 == 2]:
    violins['bodies'][i].set_facecolor('orange')

ax.set_xticks([4*x+2 for x in np.arange(indel.shape[1])])
ax.set_xticklabels(indel.columns, ha='right', rotation=30)
ax.set_ylabel('Signature Exposure Levels')
ax.set_title('INDEL Mutational Signature Exposure by ER status')
ax.legend(['Controls', 'ER Positive', 'ER Negative'], loc='upper right')

plt.tight_layout()
fig.savefig(f'/projects/b1131/saya/bbcar/plots/mutation/violin_indel_exposure_by_er_status.png')
plt.close()

