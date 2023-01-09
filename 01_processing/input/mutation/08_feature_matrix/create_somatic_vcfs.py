import numpy as np
import pandas as pd
import glob

dn = '/projects/b1131/saya/bbcar/data/02a_mutation'

sompred = pd.read_csv(f'{dn}/07_predicted_somatic/nonmatched.csv')
sompred = list(sompred.iloc[sompred.somatic.values==1,:].var_id.values)

tissue_vcfs = glob.glob(f'{dn}/02_variant_calls/tumor_only/*_DPfiltered.vcf')
tissue_sample_ids = [
    filename.split('/')[-1].split('_')[0][:-1] if filename.split('/')[-1].split('_')[0].endswith('t')
    else filename.split('/')[-1].split('_')[0] 
    for filename in tissue_vcfs
]
germline_vcfs = glob.glob(f'{dn}/02_variant_calls/germline_only/*_DPfiltered.vcf')
germline_sample_ids = [filename.split('/')[-1].split('_')[0] for filename in germline_vcfs]

tissue_only_sample_ids = set(tissue_sample_ids) - set(germline_sample_ids)

for sample_id in tissue_only_sample_ids:
    fin = f'{dn}/02_variant_calls/tumor_only/{sample_id}t_DPfiltered.vcf'
    fout = f'{dn}/07_predicted_somatic/vcfs/{sample_id}_somatic.vcf'
    # read the original calls 
    with open(fin, 'r') as f:
        lines = f.readlines()
    # write the filtered calls 
    f = open(fout, 'w')
    for line in lines:
        if line.startswith('#'):
            f.write(line)
        else:
            chrom = line.split()[0]
            pos = line.split()[1]
            ref = line.split()[3]
            alt = line.split()[4]
            var_id = f'{chrom}_{pos}_{ref}_{alt}'
            if var_id in sompred:
                f.write(line)
