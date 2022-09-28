import os
import sys
import numpy as np
import pandas as pd

din = '/projects/b1122/saya/04_cleaned_cnv'
dout = '/projects/b1122/saya/06_modified_data'

labdn = '/projects/b1122/saya' # directory where label/target is
ixdn = '/projects/b1122/saya/indices' # directory where to save train-val-test indices

#############################################################
#### First process the original, cleaned GISTIC2 outputs ####
#############################################################

for dn in os.listdir(din): # dn = {'all', 'all_ccsep', 'train_ccsep', 'train_ccsep_test_all'}
    print(f'\nWorking on directory: {dn}\n')
    # create a list of file names that we would like to reindex, based on how features are generated
    if dn=='all':
        flist=[
            'cyto_copy_conf90_all.csv',
            'gene_copy_conf90_all.csv',
            'gene_thres_conf90_all.csv',
            'gene_thres_conf90_all_plus2.csv',
            'reg_copy_conf90_all.csv',
            'reg_thres_conf90_all.csv'
        ]
    elif dn=='all_ccsep':
        flist=[
            'cyto_copy_conf90_ccsep.csv',
            'gene_copy_conf90_ccsep.csv',
            'gene_thres_conf90_ccsep.csv',
            'gene_thres_conf90_ccsep_plus2.csv',
            'reg_cc_unique_copy_gatk_all.csv'
        ]
    elif dn=='train_ccsep':
        flist=[
            'cyto_copy_conf90_ccsep.csv',
            'cyto_copy_gatk_test.csv',
            'gene_copy_conf90_ccsep.csv',
            'gene_copy_gatk_test.csv',
            'reg_cc_unique_copy_gatk_all.csv'
        ]
    elif dn=='train_ccsep_test_all':
        flist=[
            'cyto_copy_conf90_transformed.csv',
            'gene_copy_conf90_transformed.csv',
            'reg_cc_unique_copy_gatk_all.csv'
        ]
    else:
        ValueError('Unknown pattern for directory name: %s' % dn)
    # flist = os.listdir(f'{din}/{dn}')
    datadict = {}
    # Read the files into a dictionary { 'filename' : pd.DataFrame }
    for fn in flist:
        datadict[fn] = pd.read_csv(os.path.join(f'{din}/{dn}', fn), index_col=0)
    # for each file, generate a (A) study-ID indexed file and a (B) sequential integer-indexed file, save both to 'dout + dn'
    for fn in datadict.keys():
        print(f'Working on file: {fn}')
        datadict[fn]['patid'] = None
        for patid in list(datadict[fn].index):
            if type(patid)==str:
                int_id = int(patid.split('_')[0])
            else:
                int_id = patid
            datadict[fn].iloc[[x == patid for x in datadict[fn].index],-1] = int_id
        datadict[fn].set_index(['patid'], drop=True, inplace=True)
        datadict[fn].sort_index(inplace=True)
        if not os.path.isdir(f'{dout}/{dn}'):
            os.mkdir(f'{dout}/{dn}')
        datadict[fn].to_csv(os.path.join(f'{dout}/{dn}', fn.split('.')[0] + '_studyindex.csv'))
        datadict[fn].reset_index(inplace=True, drop=True) # reset in place
        datadict[fn].to_csv(os.path.join(f'{dout}/{dn}', fn.split('.')[0] + '_intindex.csv'))
        print(f'Completed reindexing for file: {fn}\n')

#################################
#### Also reindex label file ####
#################################

# load data 
labels = pd.read_csv(os.path.join(labdn, 'bbcar_label.csv'), header=0, index_col=0)

# reindex
labels['patid'] = None
for patid in list(labels.index):
    int_id = int(patid.split('_')[0])
    labels.iloc[[x == patid for x in labels.index],-1] = int_id

labels.set_index(['patid'], drop=True, inplace=True)
labels.sort_index(inplace=True)

# save with study ID as index
labels.to_csv(os.path.join(labdn, 'bbcar_label_studyindex.csv'), header=True, index=True)

# save with sequential integers as index
labels.reset_index(inplace=True, drop=True)
labels.to_csv(os.path.join(labdn, 'bbcar_label_intindex.csv'), header=True, index=True)

