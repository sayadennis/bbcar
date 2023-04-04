import glob
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
# from scipy.stats import ttest_ind
# from scipy.spatial import distance
# from scipy.cluster import hierarchy
import seaborn as sns
import matplotlib.pyplot as plt

datadir = '/projects/b1131/saya/bbcar/data'
mutdir = f'{datadir}/02a_mutation/08_feature_matrix'
# vardir = f'{datadir}/02a_mutation/02_variant_calls'
# dout_stat = '/projects/b1131/saya/bbcar/out'
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

###########################
#### Get clinical data ####
###########################

clin = pd.read_csv(f'{datadir}/00_raw/BBCaRDatabaseNU09B2_DATA_2023-04-04_0551.csv', index_col=0)

###############################
#### Calculate correlation ####
###############################

sbs1 = sig_per_sample['SBS96']['COSMIC']['SBS1']

clin.index = [int(x.split('-')[0][3:]) for x in clin.index]
clin = clin.loc[sbs1.index,:].benignbxyear

r = pearsonr(clin, sbs1)

plotdata = pd.DataFrame({'biopsy year' : clin.values, 'SBS1' : sbs1.values})

sns.lmplot(x='biopsy year', y='SBS1', data=plotdata)
plt.xlabel('Biopsy year')
plt.ylabel('SBS1 exposure')
plt.title(f'Biopsy Year & SBS1 (r={r.statistic:.2f}, p={r.pvalue:.3f})')
plt.tight_layout()
plt.savefig(f'{dout_plot}/correlation_sbs1_biopsyyear.png')
plt.close()
