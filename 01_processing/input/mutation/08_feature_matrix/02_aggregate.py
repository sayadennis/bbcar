from os.path import exists
import numpy as np
import pandas as pd
import pickle

dn = '/projects/b1131/saya/bbcar/data/02a_mutation'

data = pd.read_csv(f'{dn}/08_feature_matrix/binary_featurerows.csv', index_col=0)
anno = pd.read_csv(f'{dn}/04_ml_features/02_concat_annovar_features/annovar_features_all.csv') # annovar_features_all_bbcarpon.csv

########################################################
#### Describe mutations through annotated functions ####
########################################################

func_dictionary_fn = f'{dn}/08_feature_matrix/refGene_func_dictionary.p'

if exists(func_dictionary_fn):
    with open(func_dictionary_fn, 'rb') as f:
        func_dict = pickle.load(f)
else:
    subdata = data.iloc[data.sum(axis=1).values>1,:]
    func_dict = {}
    #### Func.refGene ####
    func_dict['Func.refGene'] = {}
    for func in anno['Func.refGene'].unique():
        func_dict['Func.refGene'][func] = []
    #
    for var_id in subdata.index:
        func = anno.iloc[anno.var_id.values==var_id,:].iloc[0]['Func.refGene']
        func_dict['Func.refGene'][func].append(var_id)
    #### ExonicFunc.refGene ####
    func_dict['ExonicFunc.refGene'] = {}
    for func in anno['ExonicFunc.refGene'].unique():
        func_dict['ExonicFunc.refGene'][func] = []
    #
    for var_id in subdata.index:
        func = anno.iloc[anno.var_id.values==var_id,:].iloc[0]['ExonicFunc.refGene']
        func_dict['ExonicFunc.refGene'][func].append(var_id)
    #### Save dictionary ####
    with open(f'{dn}/08_feature_matrix/refGene_func_dictionary.p', 'wb') as f:
        pickle.dump(func_dict, f)


# #########################
# #### Counts per gene ####
# #########################

# ct_gene = pd.DataFrame(0, index=anno['Gene.refGene'].unique(), columns=data.columns)

# for var_id in data.index:
#     gene = anno.iloc[anno.var_id.values==var_id,:]['Gene.refGene'].iloc[0]
#     var_positive_samples = data.loc[var_id].iloc[data.loc[var_id,:].values==1].index
#     ct_gene.loc[gene, var_positive_samples] += 1

# ct_gene.to_csv(f'{dn}/08_feature_matrix/cts_per_gene.csv')

#####################################################
#### Counts per gene with non-intronic mutations ####
#####################################################

ct_gene = pd.DataFrame(0, index=anno['Gene.refGene'].unique(), columns=data.columns)

for var_id in data.index:
    if var_id not in func_dict['Func.refGene']['intronic']:
        gene = anno.iloc[anno.var_id.values==var_id,:]['Gene.refGene'].iloc[0]
        var_positive_samples = data.loc[var_id].iloc[data.loc[var_id,:].values==1].index
        ct_gene.loc[gene, var_positive_samples] += 1

ct_gene.to_csv(f'{dn}/08_feature_matrix/cts_per_gene_nonintronic.csv')
