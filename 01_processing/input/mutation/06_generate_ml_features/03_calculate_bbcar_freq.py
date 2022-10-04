import sys
import glob
import numpy as np
import pandas as pd

##########################################
#### Set input and output directories ####
##########################################

din = '/projects/b1131/saya/bbcar/data/02a_mutation/04_ml_features/02_concat_annovar_features'
dout = '/projects/b1131/saya/bbcar/data/02a_mutation/04_ml_features/03_freq_added'

############################################
#### Set empty dataframe to concatenate ####
############################################

variables = [
    'var_id',
    'source',
    'sample_id',
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
    'VEST3_score',
    'CADD_raw',
    'CADD_phred',
    'GERP++_RS',
    'phyloP20way_mammalian',
    'phyloP100way_vertebrate',
    'SiPhy_29way_logOdds'
]

#### Generic PON ####

features = pd.read_csv(f'{din}/annovar_features_all.csv')

bbcar_freq = {}
for variant in features.var_id.unique():
    freq = len(features.iloc[features.var_id.values==variant,:]['sample_id'].unique())/len(features.sample_id.unique())
    bbcar_freq[variant] = freq

features['bbcar_freq'] = [bbcar_freq[x] for x in features.var_id]

features.to_csv(f'{dout}/features_annovar_bbcarfreq.csv', index=False)

#### BBCAR PON ####

features = pd.read_csv(f'{din}/annovar_features_all_bbcarpon.csv')

bbcar_freq = {}
for variant in features.var_id.unique():
    freq = len(features.iloc[features.var_id.values==variant,:]['sample_id'].unique())/len(features.sample_id.unique())
    bbcar_freq[variant] = freq

features['bbcar_freq'] = [bbcar_freq[x] for x in features.var_id]

features.to_csv(f'{dout}/features_annovar_bbcarfreq_bbcarpon.csv', index=False)
