import numpy as np
import pandas as pd
import pickle

dn = '/projects/b1131/saya/bbcar/data/02a_mutation'

data = pd.read_csv(f'{dn}/08_feature_matrix/binary_featurerows.csv', index_col=0)
anno = pd.read_csv(f'{dn}/04_ml_features/02_concat_annovar_features/annovar_features_all.csv') # annovar_features_all_bbcarpon.csv

#########################
#### Counts per gene ####
#########################

ct_gene = pd.DataFrame(0, index=anno['Gene.refGene'].unique(), columns=data.columns)

for var_id in data.index:
    gene = anno.iloc[anno.var_id.values==var_id,:]['Gene.refGene'].iloc[0]
    var_positive_samples = data.loc[var_id].iloc[data.loc[var_id,:].values==1].index
    ct_gene.loc[gene, var_positive_samples] += 1

ct_gene.to_csv(f'{dn}/08_feature_matrix/cts_per_gene.csv')

########################################################
#### Describe mutations through annotated functions ####
########################################################

subdata = data.iloc[data.sum(axis=1).values>1,:]
num = {}
num['Func.refGene'] = {}
for func in anno['Func.refGene'].unique():
    num['Func.refGene'][func] = []

for var_id in subdata.index:
    func = anno.iloc[anno.var_id.values==var_id,:].iloc[0]['Func.refGene']
    num['Func.refGene'][func].append(var_id)

num['ExonicFunc.refGene'] = {}
for func in anno['ExonicFunc.refGene'].unique():
    num['ExonicFunc.refGene'][func] = []

for var_id in subdata.index:
    func = anno.iloc[anno.var_id.values==var_id,:].iloc[0]['ExonicFunc.refGene']
    num['ExonicFunc.refGene'][func].append(var_id)

with open(f'{dn}/08_feature_matrix/refGene_func_dictionary.p', 'wb') as f:
    pickle.dump(num, f)
