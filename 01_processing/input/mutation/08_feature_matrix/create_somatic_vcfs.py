import numpy as np
import pandas as pd
import glob

dn = '/projects/b1131/saya/bbcar/data/02a_mutation'

sompred = pd.read_csv(f'{dn}/07_predicted_somatic/nonmatched.csv')
sompred = list(sompred.iloc[sompred.somatic.values==1,:].var_id.values)

tissue_only_vcf_filenames = glob.glob(f'{dn}/20221103_test/*_DPfiltered.vcf')
# tissue_only_vcf_filenames = glob.glob(f'{dn}/02_variant_calls/tumor_only/*_DPfiltered.vcf')

for filename in tissue_only_vcf_filenames:
    sample_id = filename.split('/')[-1].split('_')[0]
    fout = f'{dn}/07_predicted_somatic/vcfs/{sample_id}_somatic.vcf'
    # read the original calls 
    with open(filename, 'r') as f:
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
