import glob
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from scipy.spatial import distance
from scipy.cluster import hierarchy
import seaborn as sns
import matplotlib.pyplot as plt

datadir = '/projects/b1131/saya/bbcar/data'
mutdir = f'{datadir}/02a_mutation/08_feature_matrix'
# vardir = f'{datadir}/02a_mutation/02_variant_calls'
dout_stat = '/projects/b1131/saya/bbcar/out'
dout_plot = '/projects/b1131/saya/bbcar/plots'

######################################################
#### Get mutational signature exposure per sample ####
######################################################

dn = '/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix/signature_results'

sig_per_sample = {}

for sigtype in ['SBS96', 'DBS78', 'ID83']:
    sig_per_sample[sigtype] = {}
    samples_sig = pd.read_csv(f'{dn}/{sigtype}/Samples.txt', sep='\t', index_col=0).T
    opt_num_sig = pd.read_csv(
        f'{dn}/{sigtype}/Suggested_Solution/{sigtype}_De-Novo_Solution/Signatures/{sigtype}_De-Novo_Signatures.txt',
        sep='\t', index_col=0
    ).shape[1]
    print(f'Using {opt_num_sig} signatures for {sigtype}')
    #### De Novo Signatures ####
    opt_solution = pd.read_csv(
        f'{dn}/{sigtype}/All_Solutions/{sigtype}_{opt_num_sig}_Signatures/Signatures/{sigtype}_S{opt_num_sig}_Signatures.txt', 
        sep='\t', index_col=0
    )
    denovo_sigs_per_sample = np.dot(samples_sig, opt_solution) # realized that I had duplicate VCFs for matched samples - worknig on this now
    sig_per_sample[sigtype]['De Novo'] = pd.DataFrame(
        denovo_sigs_per_sample, 
        index=[int(x.split('_')[0]) for x in samples_sig.index], 
        columns=opt_solution.columns
    )
    #### Decomposed to COSMIC Signatures ####
    decomposed_sig = pd.read_csv(
        f'{dn}/{sigtype}/Suggested_Solution/COSMIC_{sigtype}_Decomposed_Solution/Signatures/COSMIC_{sigtype}_Signatures.txt',
        sep='\t', index_col=0
    )
    decomposed_sigs_per_sample = pd.DataFrame(
        np.dot(samples_sig, decomposed_sig),
        index=[int(x.split('_')[0]) for x in samples_sig.index],
        columns=decomposed_sig.columns
    )
    sig_per_sample[sigtype]['COSMIC'] = decomposed_sigs_per_sample

######################################
#### Get case/control assignments ####
######################################

label = pd.read_csv(f'{datadir}/clinical/bbcar_label_studyid_from_gatk_filenames.csv', index_col=0)
sample_ids = list(sig_per_sample['SBS96']['De Novo'].index)
label = label.loc[sample_ids,:]

###############################################
#### Loop through sigs and perform t-tests ####
###############################################

signatures = {}
stat_results = pd.DataFrame(columns=['Signature Type', 'Signature name', 't statistic', 'p-value'])

for sigtype in ['SBS96', 'DBS78', 'ID83']:
    #### De Novo Signatures ####
    signature = sig_per_sample[sigtype]['De Novo']
    signature_cases = signature.loc[label.iloc[label.values==1,:].index,:]
    signature_controls = signature.loc[label.iloc[label.values==0,:].index,:]
    #
    signatures[sigtype] = pd.read_csv(
        f'{mutdir}/signature_results_matched/{sigtype}/Suggested_Solution/{sigtype}_De-Novo_Solution/Signatures/{sigtype}_De-Novo_Signatures.txt', 
        sep='\t', index_col=0
    )
    #
    signature_cases_ratio = signature_cases.divide(signature_cases.sum(axis=1), axis='rows')
    signature_controls_ratio = signature_controls.divide(signature_controls.sum(axis=1), axis='rows')
    # 
    # print('\n#### Results for raw values ####')
    results = ttest_ind(signature_cases, signature_controls, axis=0)
    for signame, t, p in zip(signature.columns, results.statistic, results.pvalue):
        stat_results = pd.concat((
            stat_results, 
            pd.DataFrame({
                'Signature Type' : f'{sigtype} De Novo',
                'Signature name' : signame,
                't statistic' : t,
                'p-value' : p,
            }, index=[0])
        ))
    # 
    # print('\n#### Results for ratios ####')
    # results = ttest_ind(signature_cases_ratio, signature_controls_ratio, axis=0)
    # for signame, t, p in zip(signature.columns, results.statistic, results.pvalue):
    #     print(f'{signame}: {t:.2f} (p={p:.4f})')
    # print('\n')
    #### Decomposed to COSMIC Signatures ####
    signature = sig_per_sample[sigtype]['COSMIC']
    signature_cases = signature.loc[label.iloc[label.values==1,:].index,:]
    signature_controls = signature.loc[label.iloc[label.values==0,:].index,:]
    #
    signature_cases_ratio = signature_cases.divide(signature_cases.sum(axis=1), axis='rows')
    signature_controls_ratio = signature_controls.divide(signature_controls.sum(axis=1), axis='rows')
    # 
    # print('\n#### Results for raw values ####')
    results = ttest_ind(signature_cases, signature_controls, axis=0)
    for signame, t, p in zip(signature.columns, results.statistic, results.pvalue):
        stat_results = pd.concat((
            stat_results, 
            pd.DataFrame({
                'Signature Type' : f'{sigtype} COSMIC',
                'Signature name' : signame,
                't statistic' : t,
                'p-value' : p,
            }, index=[0])
        ))
    # 
    # print('\n#### Results for ratios ####')
    # results = ttest_ind(signature_cases_ratio, signature_controls_ratio, axis=0)
    # for signame, t, p in zip(signature.columns, results.statistic, results.pvalue):
    #     print(f'{signame}: {t:.2f} (p={p:.4f})')
    # print('\n')

stat_results.to_csv(f'{dout_stat}/20230301_matched_only_mutational_signature_exposure_ttests.csv', header=True, index=False)

###################################################################
#### Compare signature exposure levels in non-matched patients ####
###################################################################

counts = {}
sigs_per_sample_all = {}

for sigtype in ['SBS96', 'DBS78', 'ID83']:
    counts[sigtype] = pd.read_csv(f'{mutdir}/signature_results/{sigtype}/Samples.txt', sep='\t', index_col=0)
    sigs_per_sample_all[sigtype] = pd.DataFrame(
        np.dot(counts[sigtype].T, signatures[sigtype]),
        index=[int(x.split('_')[0]) for x in counts[sigtype].columns], columns=signatures[sigtype].columns
    )


label = pd.read_csv(f'{datadir}/clinical/bbcar_label_studyid_from_gatk_filenames.csv', index_col=0)

## Do statistical comparison again 
stat_results = pd.DataFrame(columns=['Signature Type', 'Signature name', 't statistic', 'p-value'])

for sigtype in ['SBS96', 'DBS78', 'ID83']:
    #### De Novo Signatures ####
    signature = sigs_per_sample_all[sigtype]
    signature_cases = signature.loc[label.iloc[label.values==1,:].index,:]
    signature_controls = signature.loc[label.iloc[label.values==0,:].index,:]
    #
    signature_cases_ratio = signature_cases.divide(signature_cases.sum(axis=1), axis='rows')
    signature_controls_ratio = signature_controls.divide(signature_controls.sum(axis=1), axis='rows')
    # 
    # print('\n#### Results for raw values ####')
    results = ttest_ind(signature_cases, signature_controls, axis=0)
    for signame, t, p in zip(signature.columns, results.statistic, results.pvalue):
        stat_results = pd.concat((
            stat_results, 
            pd.DataFrame({
                'Signature Type' : f'{sigtype} De Novo',
                'Signature name' : signame,
                't statistic' : t,
                'p-value' : p,
            }, index=[0])
        ))
