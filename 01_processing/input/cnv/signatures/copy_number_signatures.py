import glob
import itertools
import numpy as np
import pandas as pd

din = '/projects/b1131/saya/bbcar/data/02b_cnv/01_gatk_analyzed_segments'
dout = '/projects/b1131/saya/bbcar/data/02b_cnv/signatures'

###################
#### Load data ####
###################

data = pd.DataFrame(
    None,
    columns=['name', 'CONTIG', 'START', 'END', 'NUM_POINTS_COPY_RATIO', 'MEAN_LOG2_COPY_RATIO', 'CALL']
)

for fn in glob.glob(din + '/*.csv'):
    samplen = fn.split('/')[-1].split('.')[0] # example: '1004_Tissue'
    sampledata = pd.read_csv(fn) # data = GATK segment file
    sampledata = sampledata.iloc[:,:7] # only select columns necessary for GISTIC (rest is Zexian's add-ons. GISTIC throws error if kept)
    sampledata['name'] = list(itertools.repeat(samplen, len(sampledata))) # add column indicating sample name
    data = pd.concat((data, sampledata)) # append this sample to data

##########################
#### Define functions ####
##########################

def seglen_bp(line: pd.Series):
    return line['END'] - line['START']

def seg_len_category(seglen_bp):
    if ((seglen_bp > 0) & (seglen_bp <= (100 * 1e+3))):
        return '0 - 100 kb'
    elif (seglen_bp <= 1e+6):
        return '100 kb - 1 Mb'
    elif (seglen_bp <= 10 * 1e+6):
        return '1 Mb - 10 Mb'
    elif (seglen_bp <= 40 * 1e+6):
        return '10 Mb - 40 Mb'
    else:
        return '> 40 Mb'

ampdel_bins = np.quantile(2**(data.MEAN_LOG2_COPY_RATIO), [.2,.4,.6,.8,.95])

def seg_ampdel_category(mean_log2_copy_ratio, ampdel_bins=ampdel_bins):
    larger_than = np.where([mean_log2_copy_ratio > val for val in ampdel_bins])[0]
    if len(larger_than)==len(ampdel_bins):
        lower = ampdel_bins[len(ampdel_bins)-1]
        return f'> {lower:.2f}'
    elif len(larger_than)==0:
        upper = ampdel_bins[0]
        return f'<= {upper:.2f}'
    else:
        lower_ix = np.max(larger_than)
        upper_ix = lower_ix + 1
        lower = ampdel_bins[lower_ix]
        upper = ampdel_bins[upper_ix]
        return f'{lower:.2f} - {upper:.2f}'

# Add segment category by length
data['SEG CAT LENGTH'] = data.apply(seglen_bp, axis=1).apply(seg_len_category)
data['SEG CAT LEVELS'] = (2**(data.MEAN_LOG2_COPY_RATIO)).apply(seg_ampdel_category)

#########################################
#### Sample x length category matrix ####
#########################################

categories = data.groupby(['name', 'SEG CAT LENGTH', 'CALL']).size().reset_index()
counts = pd.pivot(categories, index='name', columns=('SEG CAT LENGTH', 'CALL'))

# clean up 
counts.index.name = None

tuples = list(zip(
    counts.columns.get_level_values('SEG CAT LENGTH'),
    counts.columns.get_level_values('CALL')
))
multcol = pd.MultiIndex.from_tuples(tuples, names=['SEG CAT LENGTH', 'CALL'])
counts.columns = multcol
counts = counts.iloc[:,counts.columns.get_level_values('CALL')!='0'] # Remove CALL = 0 (no amp/del)
counts.columns.names = [None, None]
counts = counts.fillna(0)

# calculate ratios 
ratios = counts.divide(counts.sum(axis=1), axis=0)
ratios['sumcounts'] = counts.sum(axis=1)

# Save 
counts.to_csv(f'{dout}/seglen_category_call_counts_per_sample.csv')
ratios.to_csv(f'{dout}/seglen_category_call_ratios_per_sample.csv')

##########################################################
#### Sample x length category x amp/del levels matrix ####
##########################################################

categories = data.groupby(['name', 'SEG CAT LENGTH', 'SEG CAT LEVELS', 'CALL']).size().reset_index()
counts = pd.pivot(categories, index='name', columns=('SEG CAT LENGTH', 'SEG CAT LEVELS', 'CALL'))

# clean up 
counts.index.name = None

tuples = list(zip(
    counts.columns.get_level_values('SEG CAT LENGTH'),
    counts.columns.get_level_values('SEG CAT LEVELS'),
    counts.columns.get_level_values('CALL')
))
multcol = pd.MultiIndex.from_tuples(tuples, names=['SEG CAT LENGTH', 'SEG CAT LEVELS', 'CALL'])
counts.columns = multcol
counts = counts.iloc[:,counts.columns.get_level_values('CALL')!='0'] # Remove CALL = 0 (no amp/del)
counts.columns.names = [None, None, None]
counts = counts.fillna(0)

# calculate ratios 
ratios = counts.divide(counts.sum(axis=1), axis=0)
ratios['sumcounts'] = counts.sum(axis=1)

# Save 
counts.to_csv(f'{dout}/seglen_ampdel_category_call_counts_per_sample.csv')
ratios.to_csv(f'{dout}/seglen_ampdel_category_call_ratios_per_sample.csv')
