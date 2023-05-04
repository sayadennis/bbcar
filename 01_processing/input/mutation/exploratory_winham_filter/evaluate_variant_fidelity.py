import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt

dout = '/projects/b1131/saya/bbcar/plots/mutation'

#########################
#### Called Variants ####
#########################

pon_source = 'bbcar'

germline_vcfs = glob.glob(f'/projects/b1131/saya/bbcar/data/02a_mutation/02_variant_calls/germline_only/*_DPfiltered.vcf')
germline_sample_ids = [filename.split('/')[-1].split('_')[0] for filename in germline_vcfs]

called = {}
for filter_type in ['liberal', 'classical', 'strict']:
    called[filter_type] = {}
    din = f'/projects/b1131/saya/bbcar/exploratory_winham_filter/{filter_type}/02_variant_calls/tumor_only'
    # dout = f'/projects/b1131/saya/bbcar/exploratory_winham_filter/{filter_type}/evaluate_calls'
    for sample_id in germline_sample_ids:
        called[filter_type][sample_id] = set()
        if len(glob.glob(f'{din}/{sample_id}_*'))==0:
            continue
        fin = glob.glob(f'{din}/{sample_id}_*')[0]
        with open(fin, 'r') as f:
            lines = f.readlines()
        #
        for line in lines:
            if not line.startswith('#'):
                chrom = line.split()[0]
                pos = line.split()[1]
                ref = line.split()[3]
                alt = line.split()[4]
                var_id = f'{chrom}_{pos}_{ref}_{alt}'
                called[filter_type][sample_id].add(var_id)

ground_truth = {}
for filter_type in ['liberal', 'classical', 'strict']:
    ground_truth[filter_type] = {}
    for sample_id in germline_sample_ids:
        if sample_id not in called['liberal'].keys():
            continue
        ground_truth[filter_type][sample_id] = set()
        din = f'/projects/b1131/saya/bbcar/exploratory_winham_filter/{filter_type}/02_variant_calls/tumor_normal'
        # dout = f'/projects/b1131/saya/bbcar/exploratory_winham_filter/{filter_type}/evaluate_calls'
        if len(glob.glob(f'{din}/{sample_id}_*'))==0:
            continue
        fin = glob.glob(f'{din}/{sample_id}_*')[0]
        with open(fin, 'r') as f:
            lines = f.readlines()
        #
        for line in lines:
            if not line.startswith('#'):
                chrom = line.split()[0]
                pos = line.split()[1]
                ref = line.split()[3]
                alt = line.split()[4]
                var_id = f'{chrom}_{pos}_{ref}_{alt}'
                ground_truth[filter_type][sample_id].add(var_id)

## There are two samples where calls are zero. Remove these samples 
zero_calls = []
for sample_id in germline_sample_ids:
    if len(called['strict'][sample_id])==0:
        zero_calls.append(sample_id)

for sample_id in zero_calls:
    germline_sample_ids.remove(sample_id)

performance = pd.DataFrame(columns=['FILTER', 'F1', 'TP', 'FP', 'FN', 'TP_STD', 'FP_STD', 'FN_STD'])
for filter_type in ['liberal', 'classical', 'strict']:
    tp = [len(called[filter_type][sample_id].intersection(ground_truth[filter_type][sample_id])) for sample_id in germline_sample_ids]
    tp_mean = np.mean(tp)
    tp_std = np.std(tp)
    fp = [len(called[filter_type][sample_id] - ground_truth[filter_type][sample_id]) for sample_id in germline_sample_ids]
    fp_mean = np.mean(fp)
    fp_std = np.std(fp)
    fn = [len(ground_truth[filter_type][sample_id] - called[filter_type][sample_id]) for sample_id in germline_sample_ids]
    fn_mean = np.mean(fn)
    fn_std = np.std(fn)
    f1 = tp_mean / (tp_mean + 1/2 * (fp_mean + fn_mean))
    performance = pd.concat((
        performance, 
        pd.DataFrame(
            [filter_type, f1, tp_mean, fp_mean, fn_mean, tp_std, fp_std, fn_std],
            columns=[0], index=['FILTER', 'F1', 'TP', 'FP', 'FN', 'TP_STD', 'FP_STD', 'FN_STD']
        ).T
    ), axis=0)

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(12,6))

## Calls that are somatic 
ax[0].violinplot(
    [[len(called[filter_type][sample_id].intersection(ground_truth[filter_type][sample_id])) for sample_id in germline_sample_ids] for filter_type in ['liberal', 'classical', 'strict']],
    positions=[0,1,2]
)
ax[0].set_xticks(np.arange(len(performance['FILTER'].unique())))
ax[0].set_xticklabels(performance['FILTER'].unique(), rotation=45, ha='right')
ax[0].set_ylabel('Number of calls')
ax[0].set_ylim(bottom=0)
# ax[0].set_ylim([0, np.max(performance['TP']*1.01)+np.max(performance['TP_STD'])]) # np.min(performance['TP'])*0.99
ax[0].set_title('Somatic variants called')

## Calls that are not somatic
ax[1].violinplot(
    [[len(called[filter_type][sample_id] - ground_truth[filter_type][sample_id]) for sample_id in germline_sample_ids] for filter_type in ['liberal', 'classical', 'strict']],
    positions=[0,1,2]
)
ax[1].set_xticks(np.arange(len(performance['FILTER'].unique())))
ax[1].set_xticklabels(performance['FILTER'].unique(), rotation=45, ha='right')
ax[1].set_ylabel('Number of calls')
ax[1].set_ylim(bottom=0)
ax[1].set_title('Non-somatic variants called')

## False negative counts
ax[2].violinplot(
    [[len(ground_truth[filter_type][sample_id] - called[filter_type][sample_id]) for sample_id in germline_sample_ids] for filter_type in ['liberal', 'classical', 'strict']],
    positions=[0,1,2]
)
ax[2].set_xticks(np.arange(len(performance['FILTER'].unique())))
ax[2].set_xticklabels(performance['FILTER'].unique(), rotation=45, ha='right')
ax[2].set_ylabel('Number of calls')
ax[2].set_ylim(bottom=0)
# ax[2].set_ylim([0, np.max(performance['FN']*1.02)]) # np.min(performance['FN'])*0.98
ax[2].set_title('Somatic variants NOT called')

fig.suptitle('Comparing Called Variants by Type of Winham Filters', fontsize=16)
plt.tight_layout()
fig.savefig(f'{dout}/violinplots_var_call_fidelity_winham_filters.png')
plt.close()

# #######################################
# #### Variants predicted as somatic ####
# #######################################

# # sompred = pd.read_csv(f'{dn}/07_predicted_somatic/nonmatched.csv')
# # sompred = list(sompred.iloc[sompred.somatic.values==1,:].var_id.values)

# pon_source = 'bbcar'

# germline_vcfs = glob.glob(f'/projects/b1131/saya/bbcar/data/02a_mutation/02_variant_calls/germline_only/*_DPfiltered.vcf')
# germline_sample_ids = [filename.split('/')[-1].split('_')[0] for filename in germline_vcfs]

# pred_som = {}
# for filter_type in ['liberal', 'classical', 'strict']:
#     pred_som[filter_type] = {}
#     din = f'/projects/b1131/saya/bbcar/exploratory_winham_filter/{filter_type}/02_variant_calls/tumor_only'
#     # dout = f'/projects/b1131/saya/bbcar/exploratory_winham_filter/{filter_type}/evaluate_calls'
#     for sample_id in germline_sample_ids:
#         pred_som[filter_type][sample_id] = set()
#         if len(glob.glob(f'{din}/{sample_id}_*'))==0:
#             continue
#         fin = glob.glob(f'{din}/{sample_id}_*')[0]
#         with open(fin, 'r') as f:
#             lines = f.readlines()
#         #
#         for line in lines:
#             if not line.startswith('#'):
#                 chrom = line.split()[0]
#                 pos = line.split()[1]
#                 ref = line.split()[3]
#                 alt = line.split()[4]
#                 var_id = f'{chrom}_{pos}_{ref}_{alt}'
#                 called[filter_type][sample_id].add(var_id)

