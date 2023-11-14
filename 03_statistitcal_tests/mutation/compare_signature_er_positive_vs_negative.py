import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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

for sigcat in ['SBS96', 'DBS78', 'ID83']:
    samples = pd.read_csv(f'{din}/{sigcat}/Samples.txt', sep='\t', index_col=0).T
    denovo = pd.read_csv(f'{din}/{sigcat}/Suggested_Solution/{sigcat}_De-Novo_Solution/Signatures/{sigcat}_De-Novo_Signatures.txt', sep='\t', index_col=0)
    denovo_exposure = pd.DataFrame(np.dot(samples.values, denovo.values), index=[int(x.split('_')[0]) for x in samples.index], columns=denovo.columns)
    data['De Novo'][sigcat] = denovo_exposure

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

sbs = data['De Novo']['SBS96']

fig, ax = plt.subplots(figsize=(12,4))

pos = [1,2,3,5,6,7,9,10,11,13,14,15,17,18,19,21,22,23]
violins = ax.violinplot(
    [sbs.iloc[labels_all['ER'].values=='Control',:]['SBS96A'].values,
    sbs.iloc[labels_all['ER'].values=='0% Negative',:]['SBS96A'].values,
    sbs.iloc[labels_all['ER'].values=='10% Positive',:]['SBS96A'].values,

    sbs.iloc[labels_all['ER'].values=='Control',:]['SBS96B'].values,
    sbs.iloc[labels_all['ER'].values=='0% Negative',:]['SBS96B'].values,
    sbs.iloc[labels_all['ER'].values=='10% Positive',:]['SBS96B'].values,

    sbs.iloc[labels_all['ER'].values=='Control',:]['SBS96C'].values,
    sbs.iloc[labels_all['ER'].values=='0% Negative',:]['SBS96C'].values,
    sbs.iloc[labels_all['ER'].values=='10% Positive',:]['SBS96C'].values,

    sbs.iloc[labels_all['ER'].values=='Control',:]['SBS96D'].values,
    sbs.iloc[labels_all['ER'].values=='0% Negative',:]['SBS96D'].values,
    sbs.iloc[labels_all['ER'].values=='10% Positive',:]['SBS96D'].values,

    sbs.iloc[labels_all['ER'].values=='Control',:]['SBS96E'].values,
    sbs.iloc[labels_all['ER'].values=='0% Negative',:]['SBS96E'].values,
    sbs.iloc[labels_all['ER'].values=='10% Positive',:]['SBS96E'].values,

    sbs.iloc[labels_all['ER'].values=='Control',:]['SBS96F'].values,
    sbs.iloc[labels_all['ER'].values=='0% Negative',:]['SBS96F'].values,
    sbs.iloc[labels_all['ER'].values=='10% Positive',:]['SBS96F'].values],
    pos, widths=0.7, showmedians=True, showextrema=True
)

for i in [0,3,6,9,12,15]:
    violins['bodies'][i].set_facecolor('blue')

for i in [1,4,7,10,13,16]:
    violins['bodies'][i].set_facecolor('green')

for i in [2,5,8,11,14,17]:
    violins['bodies'][i].set_facecolor('orange')

ax.set_xticks([2,6,10,14,18,22])
ax.set_xticklabels(['SBS96A', 'SBS96B', 'SBS96C', 'SBS96D', 'SBS96E', 'SBS96F'], ha='right', rotation=30)
ax.set_ylabel('Signature Exposure Levels')
ax.set_title('SBS Mutational Signature Exposure by ER status')
ax.legend(['Controls', 'ER Positive', 'ER Negative'], loc='upper right')

plt.tight_layout()
fig.savefig(f'/projects/b1131/saya/bbcar/plots/mutation/violin_sbs_exposure_by_er_status.png')
plt.close()

## DBS

dbs = data['De Novo']['DBS78']

fig, ax = plt.subplots(figsize=(3,4))

pos = [1,2,3,5,6,7]
violins = ax.violinplot(
    [dbs.iloc[labels_all['ER'].values=='Control',:]['DBS78A'].values,
    dbs.iloc[labels_all['ER'].values=='0% Negative',:]['DBS78A'].values,
    dbs.iloc[labels_all['ER'].values=='10% Positive',:]['DBS78A'].values,

    dbs.iloc[labels_all['ER'].values=='Control',:]['DBS78B'].values,
    dbs.iloc[labels_all['ER'].values=='0% Negative',:]['DBS78B'].values,
    dbs.iloc[labels_all['ER'].values=='10% Positive',:]['DBS78B'].values],
    pos, widths=0.7, showmedians=True, showextrema=True
)

for i in [0,3]:
    violins['bodies'][i].set_facecolor('blue')

for i in [1,4]:
    violins['bodies'][i].set_facecolor('green')

for i in [2,5]:
    violins['bodies'][i].set_facecolor('orange')

ax.set_xticks([2,6])
ax.set_xticklabels(['DBS78A', 'DBS78B'], ha='right', rotation=30)
ax.set_ylabel('Signature Exposure Levels')
ax.set_title('DBS Exposure')
#ax.legend(['Controls', 'ER Positive', 'ER Negative'], loc='upper right')

plt.tight_layout()
fig.savefig(f'/projects/b1131/saya/bbcar/plots/mutation/violin_dbs_exposure_by_er_status.png')
plt.close()

## ID83

indel = data['De Novo']['ID83']

fig, ax = plt.subplots(figsize=(10,4))

pos = [1,2,3, 5,6,7, 9,10,11, 13,14,15, 17,18,19]
violins = ax.violinplot(
    [indel.iloc[labels_all['ER'].values=='Control',:]['ID83A'].values,
    indel.iloc[labels_all['ER'].values=='0% Negative',:]['ID83A'].values,
    indel.iloc[labels_all['ER'].values=='10% Positive',:]['ID83A'].values,

    indel.iloc[labels_all['ER'].values=='Control',:]['ID83B'].values,
    indel.iloc[labels_all['ER'].values=='0% Negative',:]['ID83B'].values,
    indel.iloc[labels_all['ER'].values=='10% Positive',:]['ID83B'].values,

    indel.iloc[labels_all['ER'].values=='Control',:]['ID83C'].values,
    indel.iloc[labels_all['ER'].values=='0% Negative',:]['ID83C'].values,
    indel.iloc[labels_all['ER'].values=='10% Positive',:]['ID83C'].values,

    indel.iloc[labels_all['ER'].values=='Control',:]['ID83D'].values,
    indel.iloc[labels_all['ER'].values=='0% Negative',:]['ID83D'].values,
    indel.iloc[labels_all['ER'].values=='10% Positive',:]['ID83D'].values,

    indel.iloc[labels_all['ER'].values=='Control',:]['ID83E'].values,
    indel.iloc[labels_all['ER'].values=='0% Negative',:]['ID83E'].values,
    indel.iloc[labels_all['ER'].values=='10% Positive',:]['ID83E'].values],
    pos, widths=0.7, showmedians=True, showextrema=True
)

for i in [0,3,6,9,12]:
    violins['bodies'][i].set_facecolor('blue')

for i in [1,4,7,10,13]:
    violins['bodies'][i].set_facecolor('green')

for i in [2,5,8,11,14]:
    violins['bodies'][i].set_facecolor('orange')

ax.set_xticks([2,6,10,14,18])
ax.set_xticklabels(['ID83A', 'ID83B', 'ID83C', 'ID83D', 'ID83E'], ha='right', rotation=30)
ax.set_ylabel('Signature Exposure Levels')
ax.set_title('INDEL Mutational Signature Exposure by ER status')
ax.legend(['Controls', 'ER Positive', 'ER Negative'], loc='upper right')

plt.tight_layout()
fig.savefig(f'/projects/b1131/saya/bbcar/plots/mutation/violin_indel_exposure_by_er_status.png')
plt.close()

