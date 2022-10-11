import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_theme(color_codes=True)

din = '/projects/b1131/saya/bbcar/data/02a_mutation'
dout = '/projects/b1131/saya/bbcar/plots'

data = pd.read_csv(f'{din}/08_feature_matrix/binary_featurerows.csv', index_col=0)
anno = pd.read_csv(f'{din}/04_ml_features/02_concat_annovar_features/annovar_features_all.csv') # annovar_features_all_bbcarpon.csv
ct_gene = pd.read_csv(f'{din}/08_feature_matrix/cts_per_gene.csv', index_col=0)

#######################
#### Distributions ####
#######################

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8,4))
cts = data.sum(axis=1).values
cts = np.delete(cts, np.where(cts<=1))
ax[0].hist(cts, bins=24)
ax[0].set_xlabel('Number of samples')
ax[0].set_ylabel('Mutation counts')
ax[0].set_title('Mutation Frequency Distribution')

ax[1].hist(data.sum(axis=0).values, bins=24)
ax[1].set_xlabel('Number of mutations')
ax[1].set_ylabel('Sample counts')
ax[1].set_title('Number of somatic mutation per sample')

plt.tight_layout()
fig.savefig(f'{dout}/histograms_sample_mutation_distribution.png')
plt.close()

#################################
#### Overall cluster heatmap ####
#################################

## Remove mutations that appear in 0 or 1 samples
subdata = data.iloc[data.sum(axis=1).values>20,:]

## Plot clustering heat map 
g = sns.clustermap(subdata, xticklabels=False, yticklabels=False)
g.savefig(f'{dout}/cluster_heatmap_binary.png')
plt.close()

#########################
#### Counts per gene ####
#########################

## Remove mutations that appear in 0 or 1 samples
subdata = ct_gene.iloc[(ct_gene>0).sum(axis=1).values>=20,:]

## Plot clustering heat map 
g = sns.clustermap(subdata, xticklabels=False, yticklabels=False)
g.savefig(f'{dout}/cluster_heatmap_genewise.png')
plt.close()

########################################################
#### Describe mutations through annotated functions ####
########################################################

subdata = data.iloc[data.sum(axis=1).values>1,:]
with open(f'{din}/08_feature_matrix/refGene_func_dictionary.p', 'rb') as f:
    num = pickle.load(f)

# Plot bar graphs of exonic functions and general functions ### WORK ON THIS ### 

# Re-plot heatmap with only non-intronic mutations 
subdata = data.drop(num['Func.refGene']['intronic'], axis=0)

g = sns.clustermap(subdata, xticklabels=False, yticklabels=False)
g.savefig(f'{dout}/cluster_heatmap_binary_nonintronic.png')
plt.close()

#############################################
#### Ones with higher deleterious scores ####
#############################################

subdata = data.iloc[data.sum(axis=1).values>1,:]

scores = [
    'SIFT_score',
    'Polyphen2_HDIV_score',
    'Polyphen2_HVAR_score',
    'LRT_score',
    'MutationTaster_score',
    'MutationAssessor_score',
    'FATHMM_score',
    'MetaSVM_score',
    'MetaLR_score',
    'VEST3_score',
    'CADD_raw',
    'CADD_phred',
]

# of the mutations in subdata, how many have ANY of the above scores available? 
ct_none_avail = 0
for var_id in subdata.index:
    if np.all(pd.isnull(anno.iloc[anno.var_id.values==var_id,:][scores])):
        ct_none_avail += 1 

print(f'{ct_none_avail}/{data.shape[0]} mutations have NONE of the below scores available:')
for score_name in scores:
    print(score_name)
