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
vardir = f'{datadir}/02a_mutation/02_variant_calls'
dout_stat = '/projects/b1131/saya/bbcar/out'
dout_plot = '/projects/b1131/saya/bbcar/plots'

#######################################################
#### Identify samples that might be adding effects ####
#######################################################

with open('/projects/b1131/saya/bbcar/data/sample_ids_uchicago.txt', 'r') as f:
    uchicago = [int(x.strip()) for x in f.readlines()]

sbs96 = pd.read_csv(f'{mutdir}/20230423_signature_results/SBS96/Samples.txt', sep='\t', index_col=0)
sbs96.columns = [x.split('_')[0] for x in sbs96]

# col_linkage = hierarchy.linkage(
#     distance.pdist(sbs96.to_numpy().T), method='average')

# samples_outlier_total_vars = hierarchy.cut_tree(col_linkage, 2).ravel().astype(bool)
# samples_outlier_total_vars = list([int(x) for x in sbs96.loc[:,samples_outlier_total_vars].columns])

# exclude_samples = samples_outlier_total_vars # or can be "clustered" 

######################################
#### Get case/control assignments ####
######################################

label = pd.read_csv(f'{datadir}/clinical/bbcar_label_studyid_from_gatk_filenames.csv', index_col=0)
redcap_label = pd.read_csv(f'{datadir}/clinical/bbcar_redcap_label_studyid.csv', index_col=0)
sample_ids = pd.read_csv(f'{mutdir}/denovo_signature_per_sample_SBS96.csv', index_col=0).index
for sample_id in sample_ids:
    if sample_id not in label.index:
        label = pd.concat((
            label, 
            pd.DataFrame([redcap_label.loc[sample_id,'cc_id']], index=[sample_id], columns=['label'])
        ))

#################################################
#### Visualize the nubmer of total mutations ####
#################################################

for pon_source in ['bbcar']: # , '1000g'

    data = pd.DataFrame(index=sample_ids, columns=['total called', 'predicted somatic', 'case/control', 'matched germline'])

    #### Total called variants, prior to selecting somatic ####
    total_called_tissue_normal = {}
    for vcf in glob.glob(f'{vardir}/tumor_normal/*_DPfiltered_{pon_source}pon.vcf'):
        sample_id = int(vcf.split('/')[-1].split('_')[0])
        with open(vcf, 'r') as f:
            lines = f.readlines()
        total_called_tissue_normal[sample_id] = np.invert([line.startswith('#') for line in lines]).sum()

    total_called_tissue_only = {}
    for vcf in glob.glob(f'{vardir}/tumor_only/*_DPfiltered_{pon_source}pon.vcf'):
        sample_id = int(vcf.split('/')[-1].split('_')[0])
        with open(vcf, 'r') as f:
            lines = f.readlines()
        total_called_tissue_only[sample_id] = np.invert([line.startswith('#') for line in lines]).sum()

    # for these cases, first fill in the tissue-normal matched calls and fill in the rest with ML predictions
    for key in total_called_tissue_normal.keys():
        data.loc[key,'total called'] = total_called_tissue_normal[key]

    for key in total_called_tissue_only:
        if pd.isnull(data.loc[key,'total called']):
            data.loc[key,'total called'] = total_called_tissue_only[key]

    #### Predicted somatic variants ####
    somatic = {}
    # for vcf in glob.glob(f'{datadir}/02a_mutation/07_predicted_somatic/vcfs/*.vcf'):
    for vcf in glob.glob(f'{datadir}/02a_mutation/07_predicted_somatic/vcfs/{pon_source}/*.vcf'):
        sample_id = int(vcf.split('/')[-1].split('_')[0])
        with open(vcf, 'r') as f:
            lines = f.readlines()
        somatic[sample_id] = np.invert([line.startswith('#') for line in lines]).sum()

    for sample_id in sample_ids:
        data.loc[sample_id,'predicted somatic'] = somatic[sample_id]

    #### Fill in categorical variables about case/control and matched germline status ####
    for sample_id in data.index:
        data.loc[sample_id,'case/control'] = 'case' if label.loc[sample_id,'label'] else 'control'

    for sample_id in data.index:
        data.loc[sample_id,'matched germline'] = 'matched' if sample_id in total_called_tissue_normal.keys() else 'non-matched'

    data['matched germline'] = data['matched germline'].astype('category')
    data['case/control'] = data['case/control'].astype('category')
    data['total called'] = data['total called'].astype('float64')
    data['predicted somatic'] = data['predicted somatic'].astype('float64')

    # data = data.drop(exclude_samples, axis=0)

    ######################################################
    #### Statistical test and plotting visualizations ####
    ######################################################

    vio = sns.violinplot(
        data=data, x='matched germline', y='total called', hue='case/control', 
        split=True, cut=0, scale='count'
    )
    # vio.fig.suptitle('Number of variants called by Mutect2')
    fig = vio.get_figure()
    fig.savefig(f'{dout_plot}/violin_num_vars_compare_casecontrol_matched_total_{pon_source}.png') 
    plt.close()

    vio = sns.violinplot(
        data=data, x='matched germline', y='predicted somatic', hue='case/control', 
        split=True, cut=0, scale='count'
    )
    # vio.fig.suptitle('Number of predicted somatic variants')
    fig = vio.get_figure()
    fig.savefig(f'{dout_plot}/violin_num_vars_compare_casecontrol_matched_somatic_{pon_source}.png')
    plt.close()

    ## Test whether total counts differ significantly 
    fig, ax = plt.subplots(figsize=(3,3))
    for assignment in ['case', 'control']:
        ax.hist(
            data.iloc[data['case/control'].values==assignment,:]['total called'], 
            alpha=0.5, bins=20, range=(0,18000), label=assignment
        )

    ax.legend(loc='upper right')
    ax.set_title('Total called variants')
    fig.savefig(f'{dout_plot}/histogram_total_num_called_variants_all_{pon_source}.png')
    plt.close()

    kde = sns.kdeplot(data=data, x='total called', hue='case/control')
    fig = kde.get_figure()
    fig.savefig(f'{dout_plot}/KDE_total_num_called_variants_all_{pon_source}.png')
    plt.close()

    print(f'\n## Results for all samples combined in {pon_source}')
    print(ttest_ind(
        data.iloc[data['case/control'].values=='case',:]['total called'], 
        data.iloc[data['case/control'].values=='control',:]['total called']
    ))

    # matched samples 
    for matched_status in ['matched', 'non-matched']:
        print(f'\n## Results for {matched_status} in {pon_source} PON')
        print(ttest_ind(
            data.iloc[((data['case/control'].values=='case') & (data['matched germline'].values==matched_status)),:]['total called'], 
            data.iloc[((data['case/control'].values=='control') & (data['matched germline'].values==matched_status)),:]['total called']
        ))
        fig, ax = plt.subplots(figsize=(3,3))
        for assignment in ['case', 'control']:
            ax.hist(
                data.iloc[((data['case/control'].values==assignment) & (data['matched germline'].values==matched_status)),:]['total called'], 
                alpha=0.5, bins=20, range=(0,18000), label=assignment
            )
        #
        ax.legend(loc='upper right')
        ax.set_title(f'Called variants for {matched_status} in {pon_source} PON')
        fig.savefig(f'{dout_plot}/histogram_total_num_called_variants_{pon_source}_{matched_status}.png')
        plt.close()
        #
        kde = sns.kdeplot(data=data.iloc[data['matched germline'].values==matched_status,:], x='total called', hue='case/control')
        fig = kde.get_figure()
        fig.savefig(f'{dout_plot}/KDE_total_num_called_variants_all_{pon_source}_{matched_status}.png')
        plt.close()

    # exec(open('bbcar/repo/03_statistitcal_tests/mutation/compare_signature_case_vs_control.py').read())
