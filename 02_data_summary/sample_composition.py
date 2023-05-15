import os
import re
import glob
import numpy as np
import pandas as pd

# Load data 
label = pd.read_csv('/projects/b1122/saya/bbcar_label_studyid.csv', index_col=0)
pon = pd.read_csv('/home/srd6051/bbcar/sample_ids_pon.txt', index_col=0)

# Get labels of only PON 
pon_label = label.iloc[[i in pon.index for i in label.index],:]
# Get labels of only non-PON
nonpon_label = label.iloc[[i not in pon.index for i in label.index],:]

# Record the case/control and PON/non-PON sample composition 
composition = pd.DataFrame(0, index=['case', 'control'], columns=['tissue', 'tissue & germline'])

composition.loc['case', 'tissue & germline'] = np.sum(pon_label.values==1)
composition.loc['control', 'tissue & germline'] = np.sum(pon_label.values==0)
composition.loc['case', 'tissue'] = np.sum(nonpon_label.values==1)
composition.loc['control', 'tissue'] = np.sum(nonpon_label.values==0)

print(composition)

# Identify samples that were processed at U Chicago 

tissue = glob.glob('/projects/b1131/saya/bbcar/data/00_raw/tissue/*')
germline = glob.glob('/projects/b1131/saya/bbcar/data/00_raw/germline/*')

sequenced = {
    'chicago' : [],
    'indiana' : [],
}

for sample_dir in tissue:
    sampleid = int(sample_dir.split('/')[-1])
    fastq_fn = glob.glob(f'{sample_dir}/*.fastq.gz')
    if np.all([bool(re.search(r'_[ACTG]{7}_', x)) for x in fastq_fn]):
        sequenced['chicago'].append(sampleid)
    elif np.all(np.invert([bool(re.search(r'_[ACTG]{7}_', x)) for x in fastq_fn])):
        sequenced['indiana'].append(sampleid)
    else:
        print(f'Sample {sampleid} showing mixed results.')

# How do these compare to the list from Gannon's spreadsheet and the 35 that clustered in CN signatures? 

with open('bbcar_odd_samples.txt', 'r') as f:
    clustered = [int(x.strip()) for x in f.readlines()]

gannon = [
    1009,1074,1276,129,13,133,1389,1452,182,233,
    270,361,411,417,420,425,429,464,470,520,
    624,632,648,694,708,710,737,739,763,77,
    812,82,876,89
]
 
filenames = sequenced['chicago']

all_sampleids = list(set(clustered + gannon + filenames))

pd.DataFrame(
    [[int(x in gannon) for x in all_sampleids],
     [int(x in clustered) for x in all_sampleids],
     [int(x in filenames) for x in all_sampleids]],
    index=['All_Samples_List.xlsx','Clustered in CN sig', 'From file names'],
    columns=all_sampleids
).T.to_csv('/home/srd6051/20230309_bbcar_identifying_uchicago_samples_upsetdata.csv')
