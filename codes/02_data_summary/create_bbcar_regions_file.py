import numpy as np
import pandas as pd

fin='/projects/b1122/saya/03_gistic2_out_conf90/all_lesions.conf_90.txt'
fout='/projects/b1122/saya/04_cleaned_cnv/gistic_regions_bbcar.tsv'

data = pd.read_csv(fin, sep='\t')
widelimits = data['Wide Peak Limits'][data['Amplitude Threshold']=='Actual Copy Change Given']

ref = pd.DataFrame(None,columns=['chr', 'start', 'end', 'name', 'amp_del'])

for item in widelimits.values:
    chrn=item.split(':')[0]
    start=int(item.split('(')[0].split(':')[1].split('-')[0])
    end=int(item.split('(')[0].split(':')[1].split('-')[1])
    name=data['Unique Name'].iloc[[x==item for x in widelimits.values]].values[0]
    if 'Amplification' in name:
        amp_del = '+'
    elif 'Deletion' in name:
        amp_del = '-'
    else:
        amp_del = '.'
    ref = ref.append(pd.DataFrame([[chrn, start, end, name, amp_del]], columns=['chr', 'start', 'end', 'name', 'amp_del']))

ref.to_csv(fout, header=True, index=False, sep='\t')
