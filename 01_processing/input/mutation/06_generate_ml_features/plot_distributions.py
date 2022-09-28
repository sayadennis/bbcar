from cgi import test
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

din='/projects/b1131/saya/bbcar/04_ml_features'
dout='/projects/b1131/saya/bbcar/05_data_summary'

###################
#### Load data ####
###################

variables=[
    'var_id',
    'source',
    'sample_id',
    'Func.refGene'
    'Gene.refGene',
    'ExonicFunc.refGene',
    'Func.knownGene',
    'Gene.knownGene',
    'GeneDetail.knownGene',
    'ExonicFunc.knownGene',
    'Func.ensGene',
    'Gene.ensGene',
    'GeneDetail.ensGene',
    'ExonicFunc.ensGene',
    'AF',
    'avsnp150',
    'ExAC_ALL',
    'SIFT_score',
    'SIFT_pred',
    'Polyphen2_HDIV_score',
    'Polyphen2_HDIV_pred',
    'Polyphen2_HVAR_score',
    'Polyphen2_HVAR_pred',
    'LRT_score',
    'LRT_pred',
    'MutationTaster_score',
    'MutationTaster_pred',
    'MutationAssessor_score',
    'MutationAssessor_pred',
    'FATHMM_score',
    'FATHMM_pred',
    'MetaSVM_score',
    'MetaSVM_pred',
    'MetaLR_score',
    'MetaLR_pred',
    'VEST3_score', # some rows don't have this?? 
    'CADD_raw',
    'CADD_phred',
    'GERP++_RS',
    'phyloP20way_mammalian',
    'phyloP100way_vertebrate',
    'SiPhy_29way_logOdds'
]

categorical=[
    'var_id',
    'source',
    'sample_id',
    'Func.refGene',
    'Gene.refGene',
    'ExonicFunc.refGene',
    'Func.knownGene',
    'Gene.knownGene',
    'GeneDetail.knownGene',
    'ExonicFunc.knownGene',
    'Func.ensGene',
    'Gene.ensGene',
    'GeneDetail.ensGene',
    'ExonicFunc.ensGene',
    'avsnp150',
    'SIFT_pred',
    'Polyphen2_HDIV_pred',
    'Polyphen2_HVAR_pred',
    'LRT_pred',
    'MutationTaster_pred',
    'MutationAssessor_pred',
    'FATHMM_pred',
    'MetaSVM_pred',
    'MetaLR_pred'
]

dtype_dict={}
for colname in variables:
    if colname in categorical:
        dtype_dict[colname]=str
    else:
        dtype_dict[colname]=float

features=pd.read_csv(f'{din}/annovar_features_all.csv') # , dtype=dtype_dict

# features.AF=features.AF.replace({'0.5,0.5':0.5}).astype(float)

# #############################################################
# #### Change contents of the avsnp150 column to be binary ####
# #############################################################

# nanix=features.iloc[pd.isnull(features.avsnp150).values,:].index

# features['avsnp150'].iloc[nanix]=0
# features['avsnp150'].iloc[features.avsnp150.values!=0]=1

############################################################
#### Plot missing rate of scores based on variant types ####
############################################################

miss_rate={}
test_features=features[['Func.refGene', 'Func.knownGene', 'Func.ensGene', 'SIFT_score', 'Polyphen2_HDIV_score', 'CADD_raw']]

for func_type in ['Func.refGene', 'Func.knownGene', 'Func.ensGene']:
    unique_funcs=features[func_type].unique()
    miss_rate[func_type]=pd.DataFrame(index=unique_funcs, columns=['SIFT_score', 'Polyphen2_HDIV_score', 'CADD_raw'])
    for unique_func in unique_funcs:
        miss_rate[func_type].loc[unique_func,:]=pd.isnull(test_features.iloc[test_features[func_type].values==unique_func,:]).sum()/test_features.iloc[test_features[func_type].values==unique_func,:].shape[0]


for func_type in ['Func.refGene', 'Func.knownGene', 'Func.ensGene']:
    fig, axs = plt.subplots(ncols=3, nrows=1, figsize=(15,5))
    fig.suptitle(f'Variant annotation missing rate separated by {func_type}')
    axs[0].set_ylabel('Missing rate')
    for i in range(3):
        score=['SIFT_score', 'Polyphen2_HDIV_score', 'CADD_raw'][i]
        axs[i].bar(np.arange(miss_rate[func_type].shape[0]), miss_rate[func_type][score].values)
        axs[i].set_xticks(np.arange(len(miss_rate[func_type].index)))
        axs[i].set_xticklabels(miss_rate[func_type].index, rotation=30, ha='right', fontsize=10)
        axs[i].set_title(score)
        axs[i].set_ylim([0, 1.05])
    plt.tight_layout()
    fig.savefig(f'{dout}/missrate_by_{func_type}.png')
    plt.close()


############################################
#### Plot distributions of each feature ####
############################################

plot_features=features.drop(['source', 'sample_id', 'Gene.refGene', 'Gene.knownGene', 'GeneDetail.knownGene', 'Gene.ensGene', 'GeneDetail.ensGene'], axis=1).drop_duplicates(ignore_index=True)

fig, axs = plt.subplots(ncols=6, nrows=6, figsize=(35,20))
for i in range(1,plot_features.shape[1]):
    var_name=plot_features.columns[i]
    if var_name in categorical:
        miss_rate=100*pd.isnull(plot_features[var_name]).sum()/plot_features.shape[0]
        cts=plot_features.iloc[:,i].value_counts()
        axs[(i-1)//6,(i-1)%6].bar(np.arange(len(cts)), cts)
        axs[(i-1)//6,(i-1)%6].set_xticks(np.arange(len(cts)))
        axs[(i-1)//6,(i-1)%6].set_xticklabels(cts.index, rotation=30, ha='right') # , fontsize=14
        axs[(i-1)//6,(i-1)%6].set_title(f'{var_name} ({miss_rate:.1f}% missing)', fontsize=16)
    else:
        miss_rate=100*pd.isnull(plot_features[var_name]).sum()/plot_features.shape[0]
        axs[(i-1)//6,(i-1)%6].hist(plot_features.iloc[~pd.isnull(plot_features[var_name]).values,i].to_numpy().ravel())
        axs[(i-1)//6,(i-1)%6].set_title(f'{var_name} ({miss_rate:.1f}% missing)', fontsize=16)

plt.tight_layout()
fig.savefig(f'{dout}/feature_distributions.png')
plt.close()
