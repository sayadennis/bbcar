import os
import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

##########################################
#### Set input and output directories ####
##########################################

filter_type = sys.argv[1]

din = f'/projects/b1131/saya/bbcar/exploratory_winham_filter/{filter_type}/04_ml_features/04_imputed'
dout = f'/projects/b1131/saya/bbcar/exploratory_winham_filter/{filter_type}/04_ml_features'
dix = f'/projects/b1131/saya/bbcar/exploratory_winham_filter/{filter_type}/04_ml_features/somatic_pred_ix'

for pon_source in ['bbcar']: # , '1000g'

    ###################
    #### Load data ####
    ###################

    data = pd.read_csv(f'{din}/features_imputed_{pon_source}.csv')

    #####################################
    #### Convert to input and target ####
    #####################################

    data['somatic'] = (data.source=='tumor_normal').astype(int)

    with open(f'{dout}/tissue_normal_var_id.txt', 'w') as f:
        for var_id in list(data.iloc[data.source.values=='tumor_normal',:].var_id):
            f.write(f'{var_id}\n')

    with open(f'{dout}/tissue_only_var_id.txt', 'w') as f:
        for var_id in list(data.iloc[data.source.values=='tumor_only',:].var_id):
            f.write(f'{var_id}\n')

    with open(f'{dout}/germline_var_id.txt', 'w') as f:
        for var_id in list(data.iloc[data.source.values=='germline_only',:].var_id):
            f.write(f'{var_id}\n')


    variables=[
        'var_id',
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
        'SiPhy_29way_logOdds',
        'somatic'
    ]

    Xy_nonmatched = data.iloc[data.source.values=='tumor_only',:][variables].drop_duplicates(ignore_index=True).set_index('var_id', drop=True)
    X_nonmatched = Xy_nonmatched.iloc[:,:-1]

    Xy_matched = data.iloc[data.source.values!='tumor_only',:][variables].drop_duplicates(ignore_index=True).set_index('var_id', drop=True)
    X_matched = Xy_matched.iloc[:,:-1]
    y_matched = Xy_matched.iloc[:,-1]

    #######################################
    #### Create train and test indices ####
    #######################################

    train_ix, test_ix = train_test_split(X_matched.index, test_size=.2, random_state=43, shuffle=True)

    ##############
    #### Save ####
    ##############

    X_matched.to_csv(f'{dout}/input_matched_{pon_source}.csv', index=True, header=True)
    y_matched.to_csv(f'{dout}/target_matched_{pon_source}.csv', index=True, header=True)

    X_nonmatched.to_csv(f'{dout}/input_nonmatched_{pon_source}.csv', index=True, header=True)

    if not os.path.isdir(f'{dix}/{pon_source}'):
        os.makedirs(f'{dix}/{pon_source}')
    
    pd.DataFrame(index=train_ix).to_csv(f'{dix}/{pon_source}/train_index.txt', header=False)
    pd.DataFrame(index=test_ix).to_csv(f'{dix}/{pon_source}/test_index.txt', header=False)
