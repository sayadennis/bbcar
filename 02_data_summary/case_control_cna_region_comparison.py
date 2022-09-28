import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

din='/projects/b1122/saya'
dout='/projects/b1122/saya/05_data_summary'

cases = pd.read_csv(f'{din}/03_gistic2_out_conf90_cases/all_lesions.conf_90.txt', sep='\t')
controls = pd.read_csv(f'{din}/03_gistic2_out_conf90_controls/all_lesions.conf_90.txt', sep='\t')

## Plot a Venn diagrams for cytobands with aberrant copy numbers 

case_cyto = list(cases['Descriptor'].unique())
control_cyto = list(controls['Descriptor'].unique())

ol = np.sum([x in control_cyto for x in case_cyto])

venn2(subsets = (len(case_cyto)-ol, len(control_cyto)-ol, ol), set_labels = ('Case cytobands', 'Control cytobands'))
plt.title('Number of cytobands aberrant in cases/controls')
plt.savefig(f'{dout}/cna_cytoband_case_control_comparison.png')
plt.close()

## Evaluate region overlap 

ol_thres = 1e+3 # minimum amount of overlap between regions to consider "overlapping" 

case_regions = pd.DataFrame(
    index=cases.iloc[['CN' not in x for x in cases['Unique Name'].values],:]['Unique Name'].values, # region names 
    columns=['chr', 'start', 'end']
)
control_regions = pd.DataFrame(
    index=controls.iloc[['CN' not in x for x in controls['Unique Name'].values],:]['Unique Name'].values, # region names 
    columns=['chr', 'start', 'end']
)

for region in case_regions.index:
    peak_limits = cases.iloc[cases['Unique Name'].values==region, [x=='Wide Peak Limits' for x in cases.columns]].values[0][0].strip()
    case_regions.loc[region,'chr'] = int(peak_limits.split(':')[0].split('chr')[1]) # get the integer value of "chr<num>"
    case_regions.loc[region,'start'] = int(peak_limits.split(':')[1].split('-')[0])
    case_regions.loc[region,'end'] = int(peak_limits.split(':')[1].split('-')[1].split('(')[0])

for region in control_regions.index:
    peak_limits = controls.iloc[controls['Unique Name'].values==region, [x=='Wide Peak Limits' for x in controls.columns]].values[0][0].strip()
    control_regions.loc[region,'chr'] = int(peak_limits.split(':')[0].split('chr')[1]) # get the integer value of "chr<num>"
    control_regions.loc[region,'start'] = int(peak_limits.split(':')[1].split('-')[0])
    control_regions.loc[region,'end'] = int(peak_limits.split(':')[1].split('-')[1].split('(')[0])

ol_mx = pd.DataFrame(0, index=case_regions.index, columns=control_regions.index)
for region in case_regions.index:
    chrom = case_regions.loc[region, 'chr']
    start = case_regions.loc[region, 'start']
    end = case_regions.loc[region, 'end']
    # pick out the control_regions that are on the same chromosome for ease of handling 
    samechrom = (control_regions['chr'].values==chrom)
    start_contained = ((control_regions['start'] < start) & (control_regions['end'] > start))
    end_contained = ((control_regions['start'] < end) & (control_regions['end'] > end))
    # fill in the overlap matrix 
    # ol_mx.loc[region,:] = [np.max([control_regions['end']-start, end-control_regions['start']], axis=0)[i] if samechrom[i] else 0 for i in range(len(samechrom))]
    ol = []
    for i in range(len(samechrom)):
        if samechrom[i]:
            if (start_contained[i] & end_contained[i]):
                ol.append(end-start)
            elif (start_contained[i]):
                ol.append((control_regions['end']-start)[i]) # looking at i'th region of the control 
            elif (end_contained[i]):
                ol.append((end-control_regions['start'])[i])
            else:
                ol.append(0)
        else:
            ol.append(0)
    ol_mx.loc[region,:] = ol

olsize_control = np.sum(ol_mx.sum(axis=0).values!=0)
olsize_case = np.sum(ol_mx.sum(axis=1).values!=0)

case_unique_regions = list(ol_mx.index[ol_mx.sum(axis=1).values==0])
control_unique_regions = list(ol_mx.columns[ol_mx.sum(axis=0).values==0])

N = 2
unique = (len(case_unique_regions), len(control_unique_regions))
overlapping = (olsize_case, olsize_control)
ind = np.arange(N)  
width = 0.5

fig = plt.subplots() # figsize =(10, 7)
p1 = plt.bar(ind, unique, width)
p2 = plt.bar(ind, overlapping, width,
             bottom = unique)

plt.ylabel('Number of regions')
plt.title('Number of unique and overlapping regions for cases/controls')
plt.xticks(ind, ('Cases', 'Controls'))
# plt.yticks(np.arange(0, 81, 10))
plt.legend((p1[0], p2[0]), ('Unique', 'Overlapping'))
plt.savefig(f'{dout}/cna_region_case_control_overlap.png')
plt.close()

with open(f'{dout}/cna_region_list_case_unique.txt', 'w') as f:
    for item in case_unique_regions:
        f.write('%s\n' % item)

with open(f'{dout}/cna_region_list_control_unique.txt', 'w') as f:
    for item in control_unique_regions:
        f.write('%s\n' % item)

case_amp_genes = pd.read_csv(f'{din}/03_gistic2_out_conf90_cases/amp_genes.conf_90.txt', sep='\t', index_col=0)
case_del_genes = pd.read_csv(f'{din}/03_gistic2_out_conf90_cases/del_genes.conf_90.txt', sep='\t', index_col=0)
control_amp_genes = pd.read_csv(f'{din}/03_gistic2_out_conf90_controls/amp_genes.conf_90.txt', sep='\t', index_col=0)
control_del_genes = pd.read_csv(f'{din}/03_gistic2_out_conf90_controls/del_genes.conf_90.txt', sep='\t', index_col=0)

case_unique_amp_genes = []
case_unique_del_genes = []
control_unique_amp_genes = []
control_unique_del_genes = []

for region in case_unique_regions:
    widepeaks = 'chr%s:%s-%s' % tuple(case_regions.loc[region,:].values)
    if 'Amplification' in region:
        case_unique_amp_genes = case_unique_amp_genes + list(case_amp_genes.iloc[3:,case_amp_genes.loc['wide peak boundaries'].values==widepeaks].dropna().values.T[0])
    elif 'Deletion' in region:
        case_unique_del_genes = case_unique_del_genes + list(case_del_genes.iloc[3:,case_del_genes.loc['wide peak boundaries'].values==widepeaks].dropna().values.T[0])

for region in control_unique_regions:
    widepeaks = 'chr%s:%s-%s' % tuple(control_regions.loc[region,:].values)
    if 'Amplification' in region:
        control_unique_amp_genes = control_unique_amp_genes + list(control_amp_genes.iloc[3:,control_amp_genes.loc['wide peak boundaries'].values==widepeaks].dropna().values.T[0])
    elif 'Deletion' in region:
        control_unique_del_genes = control_unique_del_genes + list(control_del_genes.iloc[3:,control_del_genes.loc['wide peak boundaries'].values==widepeaks].dropna().values.T[0])

with open(f'{dout}/case_unique_amp_genes.txt', 'w') as f:
    for item in case_unique_amp_genes:
        f.write('%s\n' % item)

with open(f'{dout}/case_unique_del_genes.txt', 'w') as f:
    for item in case_unique_del_genes:
        f.write('%s\n' % item)

with open(f'{dout}/control_unique_amp_genes.txt', 'w') as f:
    for item in control_unique_amp_genes:
        f.write('%s\n' % item)

with open(f'{dout}/control_unique_del_genes.txt', 'w') as f:
    for item in control_unique_del_genes:
        f.write('%s\n' % item)

# venn2(subsets = (len(ol_mx.index)-ol, len(ol_mx.columns)-ol, ol), set_labels = ('Case regions', 'Control regions'))
# plt.title('Overlap Between Regions Aberrant in Cases/Controls')
# plt.savefig(f'{dout}/cna_region_case_control_comparison.png')
# plt.close()
