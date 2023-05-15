from os.path import exists
import pickle
import numpy as np
import pandas as pd
from scipy.spatial import distance
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_theme(color_codes=True)

din = '/projects/b1131/saya/bbcar/data/02a_mutation'
dout = '/projects/b1131/saya/bbcar/plots'

ct_each = pd.read_csv(f'{din}/08_feature_matrix/binary_featurerows.csv', index_col=0)
ct_gene = pd.read_csv(f'{din}/08_feature_matrix/cts_per_gene.csv', index_col=0)
ct_gene_nonintronic = pd.read_csv(f'{din}/08_feature_matrix/cts_per_gene_nonintronic.csv', index_col=0)
anno = pd.read_csv(f'{din}/04_ml_features/02_concat_annovar_features/annovar_features_all.csv') # annovar_features_all_bbcarpon.csv
target = pd.read_csv('/projects/b1131/saya/bbcar/data/02b_cnv/bbcar_label_studyid.csv', index_col=0)
encode = pd.read_csv('/projects/b1131/saya/bbcar/genome_resources/encode_regulatory/GRCh38-cCREs.bed', sep='\t', header=None)

with open(f'{din}/08_feature_matrix/refGene_func_dictionary.p', 'rb') as f:
    func_dict = pickle.load(f)

## Remove mutations that appear in 0 or 1 samples
ct_each = ct_each.iloc[ct_each.sum(axis=1).values>1,:]

############################################################################
#### Look at distributions of mutations that are in ENCODE's regulatory ####
############################################################################

regulatory_anno_fn = f'{din}/08_feature_matrix/annovar_features_regulatory.csv'
if exists(regulatory_anno_fn):
    regulatory_anno = pd.read_csv(regulatory_anno_fn)
else:
    regulatory = []
    for var_id in anno.var_id.unique():
        chrom = var_id.split('_')[0]
        pos = int(var_id.split('_')[1])
        overlap_regions = encode.iloc[(
            (encode[0].values==chrom) 
            & (encode[1].values<=pos) 
            & (encode[2].values>=pos)
        ),:]
        if overlap_regions.shape[0]>0:
            regulatory.append(var_id)
    # 
    regulatory_anno = anno.iloc[[x in regulatory for x in anno.var_id],:]
    regulatory_anno.to_csv(f'{din}/08_feature_matrix/annovar_features_regulatory.csv', index=False)

uniq_reg_vars = regulatory_anno.drop(['source', 'sample_id'], axis=1).drop_duplicates(ignore_index=True)['Func.refGene']
uniq_reg_vars.name = 'ENCODE cCRE'
uniq_vars = anno.drop(['source', 'sample_id'], axis=1).drop_duplicates(ignore_index=True)['Func.refGene']
uniq_vars.name = 'All'

pd.concat((
    uniq_reg_vars['Func.refGene'].value_counts(),
    uniq_vars['Func.refGene'].value_counts()
), axis=1)

# #######################
# #### Distributions ####
# #######################

# fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(8,4))
# cts = ct_each.sum(axis=1).values
# cts = np.delete(cts, np.where(cts<=1))
# ax[0].hist(cts, bins=24)
# ax[0].set_xlabel('Number of samples')
# ax[0].set_ylabel('Mutation counts')
# ax[0].set_title('Mutation Frequency Distribution')

# ax[1].hist(ct_each.sum(axis=0).values, bins=24)
# ax[1].set_xlabel('Number of mutations')
# ax[1].set_ylabel('Sample counts')
# ax[1].set_title('Number of somatic mutation per sample')

# plt.tight_layout()
# fig.savefig(f'{dout}/histograms_sample_mutation_distribution.png')
# plt.close()

# #################################
# #### Overall cluster heatmap ####
# #################################

# ## Plot clustering heat map 
# subdata = ct_each.iloc[ct_each.sum(axis=1).values>20,:]

# # correlations = subdata.corr()
# # correlations_array = np.asarray(subdata.corr())

# row_linkage = hierarchy.linkage(
#     distance.pdist(subdata.to_numpy()), method='average')

# col_linkage = hierarchy.linkage(
#     distance.pdist(subdata.to_numpy().T), method='average')

# temp_labels = [target.loc[int(x),'label'] if int(x) in target.index 
#                else 0 for x in subdata.columns]
# network_pal = sns.light_palette('blue', 2)
# network_lut = dict(zip([0,1], network_pal))
# network_colors = pd.Series(temp_labels).map(network_lut)

# sns.clustermap(subdata, row_linkage=row_linkage, col_linkage=col_linkage, 
#                # row_colors=network_colors, method="average",
#                col_colors=network_colors.values, 
#                figsize=(13, 13), xticklabels=False, yticklabels=False)

# plt.savefig(f'{dout}/cluster_heatmap_binary.png')
# plt.close()

#############################################################
#### Overall cluster heatmap with non-intronic mutations ####
#############################################################

## Try removing mutations that are intronic 
nonintronic = ct_each.drop(func_dict['Func.refGene']['intronic'], axis=0)

subdata = nonintronic.iloc[nonintronic.sum(axis=1).values>20,:]

row_linkage = hierarchy.linkage(
    distance.pdist(subdata.to_numpy()), method='average')

col_linkage = hierarchy.linkage(
    distance.pdist(subdata.to_numpy().T), method='average')

temp_labels = [target.loc[int(x),'label'] if int(x) in target.index 
               else 0 for x in subdata.columns]
network_pal = sns.light_palette('orange', 2)
network_lut = dict(zip([0,1], network_pal))
network_colors = pd.Series(temp_labels).map(network_lut)

sns.clustermap(subdata, row_linkage=row_linkage, col_linkage=col_linkage, 
               # row_colors=network_colors, method="average",
               col_colors=network_colors.values, 
               figsize=(13, 13), xticklabels=False, yticklabels=False)
plt.savefig(f'{dout}/cluster_heatmap_binary_nonintronic.png')
plt.close()

#########################
#### Counts per gene ####
#########################

## Remove mutations that appear in 0 or 1 samples
subdata = ct_gene.iloc[np.argsort(-1 * ct_gene.std(axis=1)).values[:2000],:]

row_linkage = hierarchy.linkage(
    distance.pdist(subdata.to_numpy()), method='average')

col_linkage = hierarchy.linkage(
    distance.pdist(subdata.to_numpy().T), method='average')

temp_labels = [target.loc[int(x),'label'] if int(x) in target.index 
               else 0 for x in subdata.columns]
network_pal = sns.light_palette('orange', 2)
network_lut = dict(zip([0,1], network_pal))
network_colors = pd.Series(temp_labels).map(network_lut)

sns.clustermap(np.log(subdata+1), row_linkage=row_linkage, col_linkage=col_linkage, 
               # row_colors=network_colors, method="average",
               col_colors=network_colors.values, 
               figsize=(13, 13), xticklabels=False, yticklabels=False)
plt.savefig(f'{dout}/cluster_heatmap_genewise.png')
plt.close()

#####################################################
#### Counts per gene with non-intronic mutations ####
#####################################################

## Remove mutations that appear in 0 or 1 samples
subdata = ct_gene_nonintronic.iloc[np.argsort(-1 * ct_gene_nonintronic.std(axis=1)).values[:2000],:]

row_linkage = hierarchy.linkage(
    distance.pdist(subdata.to_numpy()), method='average')

col_linkage = hierarchy.linkage(
    distance.pdist(subdata.to_numpy().T), method='average')

temp_labels = [target.loc[int(x),'label'] if int(x) in target.index 
               else 0 for x in subdata.columns]
network_pal = sns.light_palette('orange', 2)
network_lut = dict(zip([0,1], network_pal))
network_colors = pd.Series(temp_labels).map(network_lut)

sns.clustermap(np.log(subdata+1), row_linkage=row_linkage, col_linkage=col_linkage, 
               # row_colors=network_colors, method="average",
               col_colors=network_colors.values, 
               figsize=(13, 13), xticklabels=False, yticklabels=False)
plt.savefig(f'{dout}/cluster_heatmap_genewise_nonintronic.png')
plt.close()

# ########################################################
# #### Describe mutations through annotated functions ####
# ########################################################

# subdata = ct_each.iloc[ct_each.sum(axis=1).values>1,:]

# # Plot bar graphs of exonic functions and general functions ### WORK ON THIS ### 

# # Re-plot heatmap with only non-intronic mutations 
# subdata = ct_each.drop(func_dict['Func.refGene']['intronic'], axis=0)

# g = sns.clustermap(subdata, xticklabels=False, yticklabels=False)
# g.savefig(f'{dout}/cluster_heatmap_binary_nonintronic.png')
# plt.close()

# #############################################
# #### Ones with higher deleterious scores ####
# #############################################

# subdata = ct_each.iloc[ct_each.sum(axis=1).values>1,:]

# scores = [
#     'SIFT_score',
#     'Polyphen2_HDIV_score',
#     'Polyphen2_HVAR_score',
#     'LRT_score',
#     'MutationTaster_score',
#     'MutationAssessor_score',
#     'FATHMM_score',
#     'MetaSVM_score',
#     'MetaLR_score',
#     'VEST3_score',
#     'CADD_raw',
#     'CADD_phred',
# ]

# # of the mutations in subdata, how many have ANY of the above scores available? 
# ct_none_avail = 0
# for var_id in subdata.index:
#     if np.all(pd.isnull(anno.iloc[anno.var_id.values==var_id,:][scores])):
#         ct_none_avail += 1 

# print(f'{ct_none_avail}/{data.shape[0]} mutations have NONE of the below scores available:')
# for score_name in scores:
#     print(score_name)

###############################
#### Mutational signatures ####
###############################

din = '/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix'
dout = '/projects/b1131/saya/bbcar/plots/mutation'

denovo = pd.read_csv(f'{din}/signature_results/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt', sep='\t', index_col=0)
cosmic = pd.read_csv(f'{din}/COSMIC_v3.3.1_SBS_GRCh38.txt', sep='\t', index_col=0)
sbs96 = pd.read_csv(f'{din}/signature_results/SBS96/Samples.txt', sep='\t', index_col=0)

sbs96.columns = [x.split('_')[0] for x in sbs96]

subdata = np.log(sbs96+1)

row_linkage = hierarchy.linkage(distance.pdist(subdata.to_numpy()), method='average')
col_linkage = hierarchy.linkage(distance.pdist(subdata.to_numpy().T), method='average')

temp_labels = [target.loc[int(x),'label'] if int(x) in target.index 
               else 0 for x in subdata.columns]
network_pal = sns.light_palette('cyan', 2)
network_lut = dict(zip([0,1], network_pal))
network_colors = pd.Series(temp_labels).map(network_lut)

sns.clustermap(subdata, row_linkage=row_linkage, col_linkage=col_linkage, 
               # row_colors=network_colors, method="average",
               col_colors=network_colors.values, 
               figsize=(13, 13), xticklabels=False, yticklabels=False)
plt.savefig(f'{dout}/cluster_heatmap_sbs96.png')
plt.close()
