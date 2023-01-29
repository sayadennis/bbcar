import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

din = '/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix'
# dout = '/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix'

##########################
#### Get case/control ####
##########################

label = pd.read_csv('/projects/b1131/saya/bbcar/data/clinical/bbcar_label_studyid_from_gatk_filenames.csv', index_col=0)
redcap_label = pd.read_csv('/projects/b1131/saya/bbcar/data/clinical/bbcar_redcap_label_studyid.csv', index_col=0)
sample_ids = pd.read_csv(f'{din}/signature_per_sample_SBS96.csv', index_col=0).index
for sample_id in sample_ids:
    if sample_id not in label.index:
        label = pd.concat((
            label, 
            pd.DataFrame([redcap_label.loc[sample_id,'cc_id']], index=[sample_id], columns=['label'])
        ))

for sigtype in ['SBS96', 'DBS78', 'ID83']:
    #### De Novo Signatures ####
    print('>> De Novo Signatures <<')
    signature = pd.read_csv(f'{din}/denovo_signature_per_sample_{sigtype}.csv', index_col=0)
    signature_cases = signature.loc[label.iloc[label.values==1,:].index,:]
    signature_controls = signature.loc[label.iloc[label.values==0,:].index,:]
    #
    signature_cases_ratio = signature_cases.divide(signature_cases.sum(axis=1), axis='rows')
    signature_controls_ratio = signature_controls.divide(signature_controls.sum(axis=1), axis='rows')
    # 
    print('\n#### Results for raw values ####')
    results = ttest_ind(signature_cases, signature_controls, axis=0)
    for signame, t, p in zip(signature.columns, results.statistic, results.pvalue):
        print(f'{signame}: {t:.2f} (p={p:.4f})')
    # 
    print('\n#### Results for ratios ####')
    results = ttest_ind(signature_cases_ratio, signature_controls_ratio, axis=0)
    for signame, t, p in zip(signature.columns, results.statistic, results.pvalue):
        print(f'{signame}: {t:.2f} (p={p:.4f})')
    print('\n')
    #### Decomposed to COSMIC Signatures ####
    print('>> Signatures Decomposed to COSMIC <<')
    signature = pd.read_csv(f'{din}/cosmic_signature_per_sample_{sigtype}.csv', index_col=0)
    signature_cases = signature.loc[label.iloc[label.values==1,:].index,:]
    signature_controls = signature.loc[label.iloc[label.values==0,:].index,:]
    #
    signature_cases_ratio = signature_cases.divide(signature_cases.sum(axis=1), axis='rows')
    signature_controls_ratio = signature_controls.divide(signature_controls.sum(axis=1), axis='rows')
    # 
    print('\n#### Results for raw values ####')
    results = ttest_ind(signature_cases, signature_controls, axis=0)
    for signame, t, p in zip(signature.columns, results.statistic, results.pvalue):
        print(f'{signame}: {t:.2f} (p={p:.4f})')
    # 
    print('\n#### Results for ratios ####')
    results = ttest_ind(signature_cases_ratio, signature_controls_ratio, axis=0)
    for signame, t, p in zip(signature.columns, results.statistic, results.pvalue):
        print(f'{signame}: {t:.2f} (p={p:.4f})')
    print('\n')
