import os
import numpy as np
import pandas as pd
import glob

def gistic_features(gistic_dir, q_thres=0.05):
    # all_data_by_genes.txt file 
    gene_copy = pd.read_csv(os.path.join(gistic_dir, 'all_data_by_genes.txt'), sep='\t', index_col=0)
    gene_copy = gene_copy.T.iloc[2:,:] # remove 'Gene ID' and 'Cytoband'
    # all_thresholded.by_genes.txt file
    gene_thres = pd.read_csv(os.path.join(gistic_dir, 'all_thresholded.by_genes.txt'), sep='\t', index_col=0)
    gene_thres = gene_thres.T.iloc[2:,:] # remove 'Gene ID' and 'Cytoband'
    # all_lesions.conf_xx.txt file
    reg = pd.read_csv(glob.glob(os.path.join(gistic_dir, 'all_lesions.conf_*.txt'))[0], sep='\t', index_col=0)
    reg = reg[reg['q values'] <= q_thres] # eliminate rows with non-significant q-values
    reg_thres = reg[reg['Amplitude Threshold']!='Actual Copy Change Given']
    reg_copy = reg[reg['Amplitude Threshold']=='Actual Copy Change Given']
    reg_thres, reg_copy = reg_thres.iloc[:,8:-1], reg_copy.iloc[:,8:-1]
    return gene_copy, gene_thres, reg_copy, reg_thres

gistic_conf90 = '/projects/b1122/saya/03_gistic2_out_conf90'
dout = '/projects/b1122/saya/04_cleaned_cnv'

gene_copy, gene_thres, reg_copy, reg_thres = gistic_features(gistic_conf90)

gene_copy.to_csv(os.path.join(gistic_conf90, 'gene_copy_conf90.csv'), header=True, index=True)
gene_thres.to_csv(os.path.join(gistic_conf90, 'gene_thres_conf90.csv'), header=True, index=True)
reg_copy.to_csv(os.path.join(gistic_conf90, 'reg_copy_conf90.csv'), header=True, index=True)
reg_thres.to_csv(os.path.join(gistic_conf90, 'reg_thres_conf90.csv'), header=True, index=True)

# gistic_conf75 = '/projects/b1122/saya/gistic2_out_conf75'
# gene_copy, gene_thres, reg_copy, reg_thres = gistic_features(gistic_conf75)

# gene_copy.to_csv(os.path.join(gistic_conf75, 'gene_copy_conf75.csv'), header=True, index=True)
# gene_thres.to_csv(os.path.join(gistic_conf75, 'gene_thres_conf75.csv'), header=True, index=True)
# reg_copy.to_csv(os.path.join(gistic_conf75, 'reg_copy_conf75.csv'), header=True, index=True)
# reg_thres.to_csv(os.path.join(gistic_conf75, 'reg_thres_conf75.csv'), header=True, index=True)