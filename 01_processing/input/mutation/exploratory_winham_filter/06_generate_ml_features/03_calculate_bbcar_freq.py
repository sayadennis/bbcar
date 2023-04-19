import os
import sys
import glob
import numpy as np
import pandas as pd

##########################################
#### Set input and output directories ####
##########################################

filter_type = sys.argv[1]

din = f'/projects/b1131/saya/bbcar/exploratory_winham_filter/{filter_type}/04_ml_features/02_concat_annovar_features'
dout = f'/projects/b1131/saya/bbcar/exploratory_winham_filter/{filter_type}/04_ml_features/03_freq_added'

if not os.path.exists(dout):
    os.makedirs(dout)

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

for pon_source in ['bbcar']: # '1000g', 
    # load the concatenated annovar features 
    features = pd.read_csv(f'{din}/annovar_features_all_{pon_source}pon.csv')
    # create binary matrix where index is var_id and column is sample_id
    matrix = features[['var_id', 'sample_id']].pivot_table(index='var_id', columns='sample_id', aggfunc=lambda x: 1, fill_value=0)
    # calculate the frequency among bbcar samples and save it in a vector 
    bbcar_freq = matrix.sum(axis=1)/matrix.shape[1]
    # add this vector as a column to the existing features matrix 
    features['bbcar_freq'] = bbcar_freq.loc[features.var_id.values].values
    # save the matrix 
    features.to_csv(f'{dout}/features_annovar_bbcarfreq_{pon_source}pon.csv', index=False)
