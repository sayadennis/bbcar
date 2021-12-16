import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

dn_gatk='/projects/b1122/saya/01_gatk_analyzed_segments'
dn_gist='/projects/b1122/saya/03_gistic2_out_conf90'

## Read all GATK output files into a dictionary 
gatk = {}
for fn in glob.glob(f'{dn_gatk}/*.csv'):
    gatk[fn.split('/')[-1].split('.')[0]] = pd.read_csv(fn)

giscase = pd.read_csv(f'{dn_gist}/all_cases/all_lesions.conf_90.txt', sep='\t')
giscont = pd.read_csv(f'{dn_gist}/all_controls/all_lesions.conf_90.txt', sep='\t')

ol_thres = 1e+3 # minimum amount of overlap between regions to consider "overlapping" 

case_regions = pd.DataFrame(
    index=giscase.iloc[['CN' not in x for x in giscase['Unique Name'].values],:]['Unique Name'].values, # region names 
    columns=['chr', 'start', 'end']
)
control_regions = pd.DataFrame(
    index=giscont.iloc[['CN' not in x for x in giscont['Unique Name'].values],:]['Unique Name'].values, # region names 
    columns=['chr', 'start', 'end']
)

for region in case_regions.index:
    peak_limits = giscase.iloc[giscase['Unique Name'].values==region, [x=='Wide Peak Limits' for x in giscase.columns]].values[0][0].strip()
    case_regions.loc[region,'chr'] = int(peak_limits.split(':')[0].split('chr')[1]) # get the integer value of "chr<num>"
    case_regions.loc[region,'start'] = int(peak_limits.split(':')[1].split('-')[0])
    case_regions.loc[region,'end'] = int(peak_limits.split(':')[1].split('-')[1].split('(')[0])

for region in control_regions.index:
    peak_limits = giscont.iloc[giscont['Unique Name'].values==region, [x=='Wide Peak Limits' for x in giscont.columns]].values[0][0].strip()
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

case_unique_regions = list(ol_mx.index[ol_mx.sum(axis=1).values==0])
control_unique_regions = list(ol_mx.columns[ol_mx.sum(axis=0).values==0])

## Work on the most significant region that is unique to cases 

casen = giscase.iloc[:,9:].shape[1] # number of cases
contn = giscont.iloc[:,9:].shape[1] # number of controls
thres = 1000 # threshold for defining proximity to peak boundaries (base pairs)

print(f'Proximity threshold used: {thres}bp')

print('###############')
print('#### Cases ####')
print('###############\n')

case_record = pd.DataFrame(
    index=[f'Case {reg}' for reg in case_unique_regions],
    columns=[
    'case/control', 'amp/del', 'coordinates', 'size', '#cases w/ aberration', '#cases inclusion', '#cases start match', '#cases end match'
    ]
)

# loop through regions that were unique to cases and print output
for reg in case_unique_regions: # reg = 'Amplification Peak   4' or something like that 
    print(f'\n#### Case {reg} ####')
    case_record.loc[f'Case {reg}', 'case/control'] = 'case'
    case_record.loc[f'Case {reg}', 'amp/del'] = reg.split()[0]
    # clarify the region coordinates
    key_chr = case_regions.loc[reg,'chr']
    key_start = case_regions.loc[reg,'start']
    key_end = case_regions.loc[reg,'end']
    print(f'chr{key_chr}:{key_start}-{key_end} ({(key_end-key_start)/1000:.2f}kb)')
    case_record.loc[f'Case {reg}', 'coordinates'] = f'chr{key_chr}:{key_start}-{key_end}'
    case_record.loc[f'Case {reg}', 'size'] = f'{(key_end-key_start)/1000:.2f}kb'
    # what proportion of cases have this region amplified? 
    ab_ct = (giscase.iloc[[reg==x for x in giscase['Unique Name'].values],['Tissue' in x for x in giscase.columns]]>0).sum(axis=1).values[0]
    print(f'Number of cases with aberration in this region: {ab_ct}/{casen} ({ab_ct/casen*100:.1f}%)\n')
    case_record.loc[f'Case {reg}', '#cases w/ aberration'] = f'{ab_ct}/{casen} ({ab_ct/casen*100:.1f}%)'
    # count number of cases with matches of aberrant segments 
    ct_include = 0
    ct_startmatch = 0
    ct_endmatch = 0
    for key in gatk.keys():
        if 'Tissue' in key:
            patient_include = 0
            patient_startmatch = 0
            patient_endmatch = 0
            # see if there are segments with similar boundaries 
            for query_chr, query_start, query_end, query_call in zip(gatk[key]['CONTIG'], gatk[key]['START'], gatk[key]['END'], gatk[key]['CALL']):
                if ((query_chr==f'chr{key_chr}') & (query_start<key_start) & (query_end>key_end) & (query_call!='0')):
                    patient_include += 1
                if ((query_chr==f'chr{key_chr}') & (abs(query_start-key_start)<thres) & (query_call!='0')):
                    patient_startmatch += 1
                if ((query_chr==f'chr{key_chr}') & (abs(query_end-key_end)<thres) & (query_call!='0')):
                    patient_endmatch += 1
            if patient_include>0:
                ct_include += 1
            if patient_startmatch>0:
                ct_startmatch += 1
            if patient_endmatch>0:
                ct_endmatch += 1
    print(f'Number of cases that had segments that included this region: {ct_include} ({ct_include/casen*100:.1f}%)')
    print(f'Number of cases that had segments that had start-site proximity: {ct_startmatch} ({ct_startmatch/casen*100:.1f}%)')
    print(f'Number of cases that had segments that had end-site proximity: {ct_endmatch} ({ct_endmatch/casen*100:.1f}%)\n')
    case_record.loc[f'Case {reg}', '#cases inclusion'] = f'{ct_include} ({ct_include/casen*100:.1f}%)'
    case_record.loc[f'Case {reg}', '#cases start match'] = f'{ct_startmatch} ({ct_startmatch/casen*100:.1f}%)'
    case_record.loc[f'Case {reg}', '#cases end match'] = f'{ct_endmatch} ({ct_endmatch/casen*100:.1f}%)'

case_record.to_csv('bbcar/case_unique_region_find_breakpoints.csv', index=True, header=True)

##########

top5 = ['Amplification Peak 134', 'Amplification Peak 228', 'Amplification Peak 182', 'Deletion Peak 142', 'Amplification Peak 188']
cts = []

for i in range(len(top5)):
    # get number of patients that has significant abberration in any of peak 1,..,i 
    cts.append(np.sum(giscase.iloc[[x in top5[:i+1] for x in giscase['Unique Name'].values],9:].sum(axis=0)!=0))

plt.bar(np.arange(len(cts)), np.array(cts)/137) # 137 is the number of cases
plt.gcf().subplots_adjust(bottom=0.25)
plt.xticks(np.arange(len(cts)), labels=top5, rotation=45, ha='right', fontsize=8)
plt.ylabel('Proportion of case samples')
plt.title('Cumulative coverage of case samples with significant aberration')
plt.savefig('bbcar/region_aberration_cumulative.png')
plt.close()

##########

print('##################')
print('#### Controls ####')
print('##################\n')

control_record = pd.DataFrame(
    index=[f'Control {reg}' for reg in control_unique_regions],
    columns=[
    'case/control', 'amp/del', 'coordinates', 'size', '#controls w/ aberration', '#controls inclusion', '#controls start match', '#controls end match'
    ]
)

## Loop through regions that were unique to CONTROLS and print output
for reg in control_unique_regions: # reg = 'Amplification Peak   4' or something like that 
    print(f'\n#### Control {reg} ####')
    control_record.loc[f'Control {reg}', 'case/control'] = 'control'
    control_record.loc[f'Control {reg}', 'amp/del'] = reg.split()[0]
    # clarify the region coordinates
    key_chr = control_regions.loc[reg,'chr']
    key_start = control_regions.loc[reg,'start']
    key_end = control_regions.loc[reg,'end']
    print(f'chr{key_chr}:{key_start}-{key_end} ({(key_end-key_start)/1000:.2f}kb)')
    control_record.loc[f'Control {reg}', 'coordinates'] = f'chr{key_chr}:{key_start}-{key_end}'
    control_record.loc[f'Control {reg}', 'size'] = f'{(key_end-key_start)/1000:.2f}kb'
    # what proportion of controls have this region amplified? 
    ab_ct = (giscont.iloc[[reg==x for x in giscont['Unique Name'].values],['Control' in x for x in giscont.columns]]>0).sum(axis=1).values[0]
    print(f'Number of controls with aberration in this region: {ab_ct}/{contn} ({ab_ct/contn*100:.1f}%)\n')
    control_record.loc[f'Control {reg}', '#controls w/ aberration'] = f'{ab_ct}/{contn} ({ab_ct/contn*100:.1f}%)'
    # count number of controls with matches of aberrant segments 
    ct_include = 0
    ct_startmatch = 0
    ct_endmatch = 0
    for key in gatk.keys():
        if 'Control' in key:
            patient_include = 0
            patient_startmatch = 0
            patient_endmatch = 0
            # see if there are segments with similar boundaries 
            for query_chr, query_start, query_end, query_call in zip(gatk[key]['CONTIG'], gatk[key]['START'], gatk[key]['END'], gatk[key]['CALL']):
                if ((query_chr==f'chr{key_chr}') & (query_start<key_start) & (query_end>key_end) & (query_call!='0')):
                    patient_include += 1
                if ((query_chr==f'chr{key_chr}') & (abs(query_start-key_start)<thres) & (query_call!='0')):
                    patient_startmatch += 1
                if ((query_chr==f'chr{key_chr}') & (abs(query_end-key_end)<thres) & (query_call!='0')):
                    patient_endmatch += 1
            if patient_include>0:
                ct_include += 1
            if patient_startmatch>0:
                ct_startmatch += 1
            if patient_endmatch>0:
                ct_endmatch += 1
    print(f'Number of controls that had segments that included this region: {ct_include} ({ct_include/contn*100:.1f}%)')
    print(f'Number of controls that had segments that had start-site proximity: {ct_startmatch} ({ct_startmatch/contn*100:.1f}%)')
    print(f'Number of controls that had segments that had end-site proximity: {ct_endmatch} ({ct_endmatch/contn*100:.1f}%)\n')
    control_record.loc[f'Control {reg}', '#controls inclusion'] = f'{ct_include} ({ct_include/contn*100:.1f}%)'
    control_record.loc[f'Control {reg}', '#controls start match'] = f'{ct_startmatch} ({ct_startmatch/contn*100:.1f}%)'
    control_record.loc[f'Control {reg}', '#controls end match'] = f'{ct_endmatch} ({ct_endmatch/contn*100:.1f}%)'

control_record.to_csv('bbcar/control_unique_region_find_breakpoints.csv', index=True, header=True)
