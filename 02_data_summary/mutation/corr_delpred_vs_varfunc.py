from os.path import exists
import numpy as np
import pandas as pd
import pickle

din = '/projects/b1131/saya/bbcar/data/02a_mutation'
dout = '/projects/b1131/saya/bbcar/out'

# read data 
anno = pd.read_csv(f'{din}/04_ml_features/02_concat_annovar_features/annovar_features_all.csv')

# remove unnecessary information
anno = anno.drop(['source', 'sample_id'], axis=1)
anno = anno.drop_duplicates(ignore_index=True)

variant_functions = {
    'ExonicFunc.refGene' : ['frameshift insertion', 'frameshift deletion', 'stopgain', 'stoploss'],
    'Func.refGene' : ['exonic;splicing', 'splicing', 'ncRNA_splicing'],
}

del_preds = [
    'SIFT_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 
    'LRT_pred', 'MutationTaster_pred', 'MutationAssessor_pred', 
    'FATHMM_pred', 'MetaSVM_pred', 'MetaLR_pred',
]

for funcname in variant_functions.keys():
    for colname in del_preds:
        pd.crosstab(anno[funcname], anno[colname]).to_csv(f'{dout}/crosstab_{funcname}_{colname}.csv')
