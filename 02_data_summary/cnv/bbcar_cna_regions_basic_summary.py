import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

din='/projects/b1122/saya/03_gistic2_out_conf90'
dout='/projects/b1122/saya/05_data_summary'

amps = pd.read_csv(os.path.join(din, 'amp_genes.conf_90.txt'), sep='\t')
dels = pd.read_csv(os.path.join(din, 'del_genes.conf_90.txt'), sep='\t')

amps_bounds = amps.iloc[2,1:-1]
amps_sizes = []
for i in range(len(amps_bounds.values)):
    start = amps_bounds[i].split(':')[1].split('-')[0]
    end = amps_bounds[i].split(':')[1].split('-')[1]
    amps_sizes.append(np.absolute(int(start) - int(end)))

np.quantile(np.array(amps_sizes).ravel(), q=np.array([0, 0.25, 0.5, 0.75, 1]))

fig, ax = plt.subplots()
ax.hist(np.array(amps_sizes)/1000, bins=15)
ax.set_xlabel('Region sizes')
ax.set_ylabel('Counts')
ax.set_xticks(np.arange((2e+4)+1, step=5000))
ax.set_xticklabels(['0Mb', '5Mb', '10Mb', '15Mb', '20Mb'])
ax.set_title('Amplification Region Sizes')
fig.savefig(os.path.join(dout, 'region_summary_amp_sizes.png'))

# get sumstats for number of genes included in regions
amps_gene_num = []
for colname in list(amps.columns)[1:-1]:
    gene_num = np.sum(amps[colname].iloc[3:,].notnull())
    amps_gene_num.append(gene_num)

np.quantile(np.array(amps_gene_num), q=np.array([0, 0.25, 0.5, 0.75, 1]))

fig, ax = plt.subplots()
ax.hist(np.array(amps_gene_num), bins=15)
ax.set_xlabel('Number of genes')
ax.set_ylabel('Counts')
ax.set_title('Number of Genes in Amplification Region')
fig.savefig(os.path.join(dout, 'region_summary_amp_gene_num.png'))

## DELETIONS 

dels_bounds = dels.iloc[2,1:-1]
dels_sizes = []
for i in range(len(dels_bounds.values)):
    start = dels_bounds[i].split(':')[1].split('-')[0]
    end = dels_bounds[i].split(':')[1].split('-')[1]
    dels_sizes.append(np.absolute(int(start) - int(end)))

np.quantile(np.array(dels_sizes).ravel(), q=np.array([0, 0.25, 0.5, 0.75, 1]))

fig, ax = plt.subplots()
ax.hist(np.array(dels_sizes)/1000, bins=15)
ax.set_xlabel('Region sizes')
ax.set_ylabel('Counts')
ax.set_xticks(np.arange((2e+4)+1, step=5000))
ax.set_xticklabels(['0Mb', '5Mb', '10Mb', '15Mb', '20Mb'])
ax.set_title('Deletion Region Sizes')
fig.savefig(os.path.join(dout, 'region_summary_del_sizes.png'))

# get sumstats for number of genes included in regions
dels_gene_num = []
for colname in list(dels.columns)[1:-1]:
    gene_num = np.sum(dels[colname].iloc[3:,].notnull())
    dels_gene_num.append(gene_num)

np.quantile(np.array(dels_gene_num), q=np.array([0, 0.25, 0.5, 0.75, 1]))

fig, ax = plt.subplots()
ax.hist(np.array(dels_gene_num), bins=15)
ax.set_xlabel('Number of genes')
ax.set_ylabel('Counts')
ax.set_title('Number of Genes in Deletion Region')
fig.savefig(os.path.join(dout, 'region_summary_del_gene_num.png'))
