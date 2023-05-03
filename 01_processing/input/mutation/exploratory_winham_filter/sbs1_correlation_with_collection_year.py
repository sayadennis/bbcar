import glob
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from scipy.stats import linregress
from scipy import stats
# from scipy.stats import ttest_ind
# from scipy.spatial import distance
# from scipy.cluster import hierarchy
import seaborn as sns
import matplotlib.pyplot as plt

datadir = '/projects/b1131/saya/bbcar/data'
mutdir = f'{datadir}/02a_mutation/08_feature_matrix'
# vardir = f'{datadir}/02a_mutation/02_variant_calls'
# dout_stat = '/projects/b1131/saya/bbcar/out'
dout_plot = '/projects/b1131/saya/bbcar/plots/mutation'

######################################################
#### Get mutational signature exposure per sample ####
######################################################

sig_per_sample = {}

for filter_type in ['liberal', 'classical', 'strict']:
    print(f'#### Filter Type: {filter_type} ####')
    sig_per_sample[filter_type] = {}
    dn = f'/projects/b1131/saya/bbcar/exploratory_winham_filter/{filter_type}/08_feature_matrix/signature_results'
    for sigtype in ['SBS96', 'DBS78', 'ID83']:
        sig_per_sample[filter_type][sigtype] = {}
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
        sig_per_sample[filter_type][sigtype]['De Novo'] = pd.DataFrame(
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
        sig_per_sample[filter_type][sigtype]['COSMIC'] = decomposed_sigs_per_sample

# Also definte no-filter data 
filter_type = 'no filter'
sig_per_sample[filter_type] = {}
dn = f'/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix/signature_results'

for sigtype in ['SBS96', 'DBS78', 'ID83']:
    sig_per_sample[filter_type][sigtype] = {}
    samples_sig = pd.read_csv(f'{dn}/{sigtype}/Samples.txt', sep='\t', index_col=0).T
    opt_num_sig = pd.read_csv(
        f'{dn}/{sigtype}/Suggested_Solution/{sigtype}_De-Novo_Solution/Signatures/{sigtype}_De-Novo_Signatures.txt',
        sep='\t', index_col=0
    ).shape[1]
    #### De Novo Signatures ####
    opt_solution = pd.read_csv(
        f'{dn}/{sigtype}/All_Solutions/{sigtype}_{opt_num_sig}_Signatures/Signatures/{sigtype}_S{opt_num_sig}_Signatures.txt', 
        sep='\t', index_col=0
    )
    denovo_sigs_per_sample = np.dot(samples_sig, opt_solution) # realized that I had duplicate VCFs for matched samples - worknig on this now
    sig_per_sample[filter_type][sigtype]['De Novo'] = pd.DataFrame(
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
    sig_per_sample[filter_type][sigtype]['COSMIC'] = decomposed_sigs_per_sample

###########################
#### Get clinical data ####
###########################

clin = pd.read_csv(f'{datadir}/00_raw/BBCaRDatabaseNU09B2_DATA_2023-04-04_0551.csv', index_col=0)
clin.index = [int(x.split('-')[0][3:]) for x in clin.index]

###############################
#### Calculate correlation ####
###############################

fig, axs = plt.subplots(ncols=4, nrows=1, figsize=(14,5))

scatter_ylim = (-15, np.max(sig_per_sample['no filter']['SBS96']['COSMIC']['SBS1'])+30)

for i, filter_type in enumerate(['no filter', 'liberal', 'classical', 'strict']):
    # SBS1 levels 
    sbs1 = sig_per_sample[filter_type]['SBS96']['COSMIC']['SBS1']
    # subset clinical data
    subclin = clin.loc[sbs1.index,:].benignbxyear
    # calculate pearson correlation 
    r = pearsonr(subclin, sbs1)
    # plot 
    plotdata = pd.DataFrame({'biopsy year' : subclin.values, 'SBS1' : sbs1.values})
    axs[i].scatter(plotdata['biopsy year'], plotdata['SBS1'])
    # Fit a linear regression line to the data
    slope, intercept = np.polyfit(plotdata['biopsy year'], plotdata['SBS1'], deg=1)
    x_min, x_max = axs[i].get_xlim()
    line_x = np.array([x_min, x_max])
    line_y = slope * line_x + intercept
    axs[i].plot(line_x, line_y)
    axs[i].set_ylim(scatter_ylim)
    axs[i].set_xlabel('Biopsy year')
    if i==0:
        axs[i].set_ylabel('SBS1 exposure')
    axs[i].set_title(f'Filter: {filter_type.upper()} (r={r.statistic:.2f}, p={r.pvalue:.3f})')

fig.suptitle('Correlations of Biopsy Year & SBS1 for 3 Winham Filters')
fig.tight_layout()
fig.savefig(f'{dout_plot}/correlation_sbs1_biopsyyear_winhamfilters.png')
plt.close()
