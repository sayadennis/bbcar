import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob

dgis='/projects/b1122/saya/02_gistic2_input' # gistic input directory 
dint='/projects/b1122/gannon/bbcar/RAW_data/int_lst' # interval directory
dout='/projects/b1122/saya/05_data_summary'

## Read the "GISTIC input" i.e. all GATK segments concatenated into one table
gisin = pd.read_csv(f'{dgis}/gistic_input_all.tsv', sep='\t', header=None)

# separate by chromosome
chroms = {}
for chrom in gisin[1].unique():
    chroms[chrom] = gisin.iloc[gisin[1].values==chrom,:]

# record the possible breakpoints
possible_breakpoints = pd.DataFrame(columns=['chrom', 'pos', 'start/end', '#samples'])

for chrom in chroms.keys():
    # find coordinates that appear as segment boundaries in >=30 samples (start and end separately)
    freq_startpos = chroms[chrom][2].value_counts().iloc[chroms[chrom][2].value_counts().values>=30]
    freq_endpos = chroms[chrom][3].value_counts().iloc[chroms[chrom][3].value_counts().values>=30]
    # record these positions and number of samples onto the dataframe
    possible_breakpoints = possible_breakpoints.append(
        pd.DataFrame({
            'chrom' : chrom, 
            'pos' : freq_startpos.index, 
            'start/end' : 'start', 
            '#samples' : freq_startpos.values
        }),
    ignore_index=True)
    possible_breakpoints = possible_breakpoints.append(
        pd.DataFrame({
            'chrom' : chrom, 
            'pos' : freq_endpos.index, 
            'start/end' : 'end', 
            '#samples' : freq_endpos.values
        }),
    ignore_index=True)

# save 
possible_breakpoints.to_csv(f'{dout}/test_find_breakpoints_gatk.csv', header=True, index=True)


thres = 100 # minimum basepair distance for proximity 

## Check whether any of these boundaries match with the sequencing intervals

# first save possible breakpoints by chromosome 
chroms = {}
for chrom in possible_breakpoints['chrom'].unique():
    chroms[chrom] = possible_breakpoints.iloc[possible_breakpoints['chrom'].values==chrom,:]

v5int = pd.read_csv(
    f'{dint}/SureSelect_v5/hg38/hg38.preprocessed.interval_list', 
    sep='\t', comment='@', header=None,
    names=['chrom', 'start', 'end', 'strand', 'other']
)
v6int = pd.read_csv(
    f'{dint}/SureSelect_v6/hg38.preprocessed.interval_list', 
    sep='\t', comment='@', header=None,
    names=['chrom', 'start', 'end', 'strand', 'other']
)

int_match = pd.DataFrame(columns=['chrom', 'pos', 'start/end of segment', 'v5/v6', 'start/end of intervals'])
non_match = pd.DataFrame(columns=['chrom', 'pos', 'start/end', '#samples'])

for chrom in chroms.keys():
    v5chrom_ints = v5int.iloc[v5int['chrom'].values==chrom,:]
    # if any of the possible_breakpoints['pos'] matches 
    for i in range(chroms[chrom].shape[0]):
        start_match = v5chrom_ints['start'].values==chroms[chrom]['pos'].iloc[i]
        end_match = v5chrom_ints['end'].values==chroms[chrom]['pos'].iloc[i]
        if np.any((start_match | end_match)):
            for matchloc in np.where((start_match | end_match))[0]:
                int_match = int_match.append({
                    'chrom' : chrom,
                    'pos' : chroms[chrom]['pos'].iloc[i],
                    'start/end of segment' : chroms[chrom]['start/end'].iloc[i],
                    'v5/v6' : 'v5',
                    'start/end of intervals' : ('start' if start_match[matchloc] else 'end')
                }, ignore_index=True)

for chrom in chroms.keys():
    v6chrom_ints = v6int.iloc[v6int['chrom'].values==chrom,:]
    # if any of the possible_breakpoints['pos'] matches 
    for i in range(chroms[chrom].shape[0]):
        start_match = v6chrom_ints['start'].values==chroms[chrom]['pos'].iloc[i]
        end_match = v6chrom_ints['end'].values==chroms[chrom]['pos'].iloc[i]
        if np.any((start_match | end_match)):
            for matchloc in np.where((start_match | end_match))[0]:
                int_match = int_match.append({
                    'chrom' : chrom,
                    'pos' : chroms[chrom]['pos'].iloc[i],
                    'start/end of segment' : chroms[chrom]['start/end'].iloc[i],
                    'v5/v6' : 'v6',
                    'start/end of intervals' : ('start' if start_match[matchloc] else 'end')
                }, ignore_index=True)
        else:
            non_match = non_match.append(chroms[chrom].iloc[i,:], ignore_index=True)


int_match.to_csv(f'{dout}/test_find_interval_match_with_breakpoints.csv', header=True, index=True)
non_match.to_csv(f'{dout}/test_find_interval_nonmatch_with_breakpoints.csv', header=True, index=True)
