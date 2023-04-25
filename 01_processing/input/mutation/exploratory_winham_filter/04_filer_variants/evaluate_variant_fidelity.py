import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt

pon_source = 'bbcar'

# sompred = pd.read_csv(f'{dn}/07_predicted_somatic/nonmatched.csv')
# sompred = list(sompred.iloc[sompred.somatic.values==1,:].var_id.values)

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

performance = pd.DataFrame(columns=['FILTER', 'F1', 'TP', 'FP', 'FN'])
for filter_type in ['liberal', 'classical', 'strict']:
    tp = np.sum([len(called[filter_type][sample_id].intersection(ground_truth[filter_type][sample_id])) for sample_id in germline_sample_ids])
    fp = np.sum([len(called[filter_type][sample_id] - ground_truth[filter_type][sample_id]) for sample_id in germline_sample_ids])
    fn = np.sum([len(ground_truth[filter_type][sample_id] - called[filter_type][sample_id]) for sample_id in germline_sample_ids])
    f1 = tp / (tp + 1/2 * (fp + fn))
    performance = pd.concat((
        performance, 
        pd.DataFrame(
            [filter_type, f1, tp, fp, fn],
            columns=[0], index=['FILTER', 'F1', 'TP', 'FP', 'FN']
        ).T
    ), axis=0)

fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(12,6))

# ## F1 score 
# ax[0].bar(
#     x=np.arange(len(performance['FILTER'].unique())),
#     height=[np.mean(performance['F1'].iloc[performance['FILTER'].values==x]) for x in ['liberal', 'classical', 'strict']],
#     yerr=[np.std(performance['F1'].iloc[performance['FILTER'].values==x]) for x in ['liberal', 'classical', 'strict']],
# )
# ax[0].set_xticks(np.arange(len(performance['FILTER'].unique())))
# ax[0].set_xticklabels([f'{int(x)} samples' for x in performance['FILTER'].unique()], rotation=45, ha='right')
# ax[0].set_ylabel('F1 score')
# ax[0].set_title('F1 Score')

## Calls that are somatic 
ax[0].bar(
    x=np.arange(len(performance['FILTER'].unique())),
    height=[np.mean(performance['TP'].iloc[performance['FILTER'].values==x]) for x in ['liberal', 'classical', 'strict']],
    # yerr=[np.std(performance['TP'].iloc[performance['FILTER'].values==x]) for x in ['liberal', 'classical', 'strict']],
)
ax[0].set_xticks(np.arange(len(performance['FILTER'].unique())))
ax[0].set_xticklabels(performance['FILTER'].unique(), rotation=45, ha='right')
ax[0].set_ylabel('Number of calls')
ax[0].set_ylim([0, np.max(performance['TP']*1.01)]) # np.min(performance['TP'])*0.99
ax[0].set_title('Somatic variants called')

## Calls that are not somatic
ax[1].bar(
    x=np.arange(len(performance['FILTER'].unique())),
    height=[np.mean(performance['FP'].iloc[performance['FILTER'].values==x]) for x in ['liberal', 'classical', 'strict']],
    # yerr=[np.std(performance['FP'].iloc[performance['FILTER'].values==x]) for x in ['liberal', 'classical', 'strict']],
)
ax[1].set_xticks(np.arange(len(performance['FILTER'].unique())))
ax[1].set_xticklabels(performance['FILTER'].unique(), rotation=45, ha='right')
ax[1].set_ylabel('Number of calls')
ax[1].set_title('Non-somatic variants called')

## False negative counts
ax[2].bar(
    x=np.arange(len(performance['FILTER'].unique())),
    height=[np.mean(performance['FN'].iloc[performance['FILTER'].values==x]) for x in ['liberal', 'classical', 'strict']],
    # yerr=[np.std(performance['FN'].iloc[performance['FILTER'].values==x]) for x in ['liberal', 'classical', 'strict']],
)
ax[2].set_xticks(np.arange(len(performance['FILTER'].unique())))
ax[2].set_xticklabels(performance['FILTER'].unique(), rotation=45, ha='right')
ax[2].set_ylabel('Number of calls')
ax[2].set_ylim([0, np.max(performance['FN']*1.02)]) # np.min(performance['FN'])*0.98
ax[2].set_title('Somatic variants NOT called')

fig.suptitle('Comparing Called Variants by Type of Winham Filters', fontsize=16)
plt.tight_layout()
fig.savefig(f'prelim_performance.png')
plt.close()

# performance.to_csv(f'{dn}/')
