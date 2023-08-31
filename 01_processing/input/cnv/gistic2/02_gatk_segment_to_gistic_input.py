import os
import numpy as np
import pandas as pd
import itertools
import glob

din = '/projects/b1131/saya/bbcar/data/02b_cnv/01_gatk_analyzed_segments'
dout = '/projects/b1131/saya/bbcar/data/02b_cnv/02_gistic2_input'
dix = '/projects/b1131/saya/bbcar/data/02b_cnv/indices'

###############################
#### All samplies combined ####
###############################

allmx = pd.DataFrame(None, columns=['name', 'CONTIG', 'START', 'END', 'NUM_POINTS_COPY_RATIO', 'MEAN_LOG2_COPY_RATIO'])

for fn in glob.glob(din + '/*.csv'):
    samplen = fn.split('/')[-1].split('.')[0] # example: '1004_Tissue'
    data = pd.read_csv(fn) # data = GATK segment file
    newdata = data.iloc[:,:6] # only select columns necessary for GISTIC (rest is Zexian's add-ons. GISTIC throws error if kept)
    newdata['name'] = list(itertools.repeat(samplen, len(newdata))) # add column indicating sample name
    allmx = pd.concat((allmx, newdata)) # append this sample to allmx

allmx['MEAN_LOG2_COPY_RATIO'] = allmx['MEAN_LOG2_COPY_RATIO'] - 1 # GISTIC and GATK uses different log2 standards
allmx.to_csv(os.path.join(dout, 'gistic_input_all.tsv'), sep='\t', header=False, index=False)

###############################################
#### All cases and all controls separately ####
###############################################

casemx = pd.DataFrame(None, columns=['name', 'CONTIG', 'START', 'END', 'NUM_POINTS_COPY_RATIO', 'MEAN_LOG2_COPY_RATIO'])
controlmx = pd.DataFrame(None, columns=['name', 'CONTIG', 'START', 'END', 'NUM_POINTS_COPY_RATIO', 'MEAN_LOG2_COPY_RATIO'])

for fn in glob.glob(din + '/*.csv'):
    samplen = fn.split('/')[-1].split('.')[0]
    data = pd.read_csv(fn)
    newdata = data.iloc[:,:6] # only select the necessary columns (first 7)
    newdata['name'] = list(itertools.repeat(samplen, len(newdata)))
    if 'Tissue' in samplen:
        casemx = pd.concat((casemx, newdata))
    elif 'Control' in samplen:
        controlmx = pd.concat((controlmx, newdata))
    # combmx = pd.concat((combmx, newdata))

casemx['MEAN_LOG2_COPY_RATIO'] = casemx['MEAN_LOG2_COPY_RATIO'] - 1
controlmx['MEAN_LOG2_COPY_RATIO'] = controlmx['MEAN_LOG2_COPY_RATIO'] - 1

casemx.to_csv(os.path.join(dout, 'gistic_input_all_cases.tsv'), sep='\t', header=False, index=False)
controlmx.to_csv(os.path.join(dout, 'gistic_input_all_controls.tsv'), sep='\t', header=False, index=False)

#########################################################
#### Training cases and training controls separately ####
#########################################################

traincasemx = pd.DataFrame(None, columns=['name', 'CONTIG', 'START', 'END', 'NUM_POINTS_COPY_RATIO', 'MEAN_LOG2_COPY_RATIO'])
traincontrolmx = pd.DataFrame(None, columns=['name', 'CONTIG', 'START', 'END', 'NUM_POINTS_COPY_RATIO', 'MEAN_LOG2_COPY_RATIO'])

all_ix = pd.read_csv(f'{dix}/index_dictionary.csv')
train_ix = list(pd.read_csv(f'{dix}/train_studyid.csv', header=None)[0])

for fn in glob.glob(din + '/*.csv'):
    samplen = fn.split('/')[-1].split('.')[0]
    data = pd.read_csv(fn)
    newdata = data.iloc[:,:6] # only select the necessary columns (first 7)
    newdata['name'] = list(itertools.repeat(samplen, len(newdata)))
    if int(samplen.split('_')[0]) in train_ix:
        if 'Tissue' in samplen:
            traincasemx = pd.concat((traincasemx, newdata))
        else: # if 'Control' in samplen:
            traincontrolmx = pd.concat((traincontrolmx, newdata))

traincasemx['MEAN_LOG2_COPY_RATIO'] = traincasemx['MEAN_LOG2_COPY_RATIO'] - 1
traincontrolmx['MEAN_LOG2_COPY_RATIO'] = traincontrolmx['MEAN_LOG2_COPY_RATIO'] - 1

traincasemx.to_csv(os.path.join(dout, 'gistic_input_train_cases.tsv'), sep='\t', header=False, index=False)
traincontrolmx.to_csv(os.path.join(dout, 'gistic_input_train_controls.tsv'), sep='\t', header=False, index=False)

###################################################################################
#### (Training cases + test all) and (training controls + test all) separately ####
###################################################################################

traincasetestallmx = pd.DataFrame(None, columns=['name', 'CONTIG', 'START', 'END', 'NUM_POINTS_COPY_RATIO', 'MEAN_LOG2_COPY_RATIO'])
traincontroltestallmx = pd.DataFrame(None, columns=['name', 'CONTIG', 'START', 'END', 'NUM_POINTS_COPY_RATIO', 'MEAN_LOG2_COPY_RATIO'])

all_ix = pd.read_csv(f'{dix}/index_dictionary.csv')
train_ix = list(pd.read_csv(f'{dix}/train_studyid.csv', header=None)[0])

for fn in glob.glob(din + '/*.csv'):
    samplen = fn.split('/')[-1].split('.')[0]
    data = pd.read_csv(fn)
    newdata = data.iloc[:,:6]
    newdata['name'] = list(itertools.repeat(samplen, len(newdata)))
    if int(samplen.split('_')[0]) in train_ix:
        if 'Tissue' in samplen:
            traincasetestallmx = pd.concat((traincasetestallmx, newdata))
        elif 'Control' in samplen:
            traincontroltestallmx = pd.concat((traincontroltestallmx, newdata))
    else: # if not in training set (i.e. in test set), append to both case/control matrices
        traincasetestallmx = pd.concat((traincasetestallmx, newdata))
        traincontroltestallmx = pd.concat((traincontroltestallmx, newdata))

traincasetestallmx['MEAN_LOG2_COPY_RATIO'] = traincasetestallmx['MEAN_LOG2_COPY_RATIO'] - 1
traincontroltestallmx['MEAN_LOG2_COPY_RATIO'] = traincontroltestallmx['MEAN_LOG2_COPY_RATIO'] - 1

traincasetestallmx.to_csv(os.path.join(dout, 'gistic_input_train_cases_test_all.tsv'), sep='\t', header=False, index=False)
traincontroltestallmx.to_csv(os.path.join(dout, 'gistic_input_train_controls_test_all.tsv'), sep='\t', header=False, index=False)
