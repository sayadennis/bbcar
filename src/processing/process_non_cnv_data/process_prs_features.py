import os
import sys
import glob
import numpy as np
import pandas as pd

din = '/projects/b1122/Zexian/Alignment/BBCAR_NEW/WES_Analysis/Haplotype/VCF'
dout = '/projects/b1122/saya/additional_features'

snpdata = pd.read_excel(os.path.join(dout, 'prs_snp_mavaddat2019_cleaned.xlsx'), header=0)
sample_ids = list(pd.read_csv('bbcar/all_samples.txt', index_col=0, header=None).index)

snpfeatures = pd.DataFrame(None, index=sample_ids, columns=snpdata['SNP'].values)

for sample_id in sample_ids:
    fn = os.path.join(din, '%s_raw_snps.vcf' % sample_id)
    snplist = []
    vcf = pd.read_csv(fn, comment='#', sep='\t', header=None)
    for i in range(snpdata.shape[0]):
        chrname = 'chr' + str(snpdata['Chromosome'].iloc[i])
        pos = snpdata['Position'].iloc[i]
        alt = snpdata['Effect Allele'].iloc[i]
        if np.sum(np.array(vcf.iloc[:,0] == chrname) & np.array(vcf.iloc[:,1] == pos) & np.array(vcf.iloc[:,4] == alt)):
            snplist.append(1)
        else:
            snplist.append(0)
    snpfeatures.loc[sample_id,:] = snplist

# remove columns that are all zero 
snpfeatures = snpfeatures.iloc[:,(snpfeatures.sum(axis=0)!=0).values]

#### Save necessary information in ML-appropriate format #### 
snpfeatures.to_csv(os.path.join(dout, 'bbcar_prs_studyindex.csv'), header=True, index=True)

# save CSV with integer index too
snpfeatures.reset_index(drop=True).to_csv(os.path.join(dout, 'bbcar_prs_intindex.csv'), header=True, index=True)
