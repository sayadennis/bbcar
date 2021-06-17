import os
import sys
import numpy as np
import pandas as pd
from scipy.special import softmax
import glob
import torch
from random_permutation import get_pval

datadir = '/projects/b1122/saya/scanmap_data'
ixdir = '/projects/b1122/saya/scanmap_data'
resultsdir = 'bbcar/scanmap/results'
train_record_dir = '/projects/b1042/ClareLab/saya/train_record'

test_indices = list(pd.read_csv(os.path.join(ixdir, 'test_indices_0.1val_0.2te.csv'), header=None, index_col=0).index)
y_test = pd.read_csv(os.path.join(datadir, 'bbcar_label_intindex.csv'), header=0, index_col=0).loc[test_indices,:].to_numpy().ravel()

##################################################
#### First save predictions of ScanMap models ####
##################################################

for fn in glob.glob(resultsdir, '/results_*.txt'):
    # find best performing model's parameters
    try:
        results = pd.read_csv(fn, skiprows=2)
        bestpf = results.sort_values('val bal acc', ascending=False).iloc[0,:]
        nc, C, wcls, lr = int(bestpf['nc']), bestpf['C'], bestpf['wcls'], bestpf['lr']
        if wcls >= 1: 
            wcls = int(wcls)
        if C >=1:
            C = int(C)
        # get feature names 
        gen_feature = fn.split('/')[-1].split('_')[1]
        cf_name = fn[fn.index('_', 30)+1:].split('.')[0] # this will be something like 'clin_prs_driversomatic'
        cf_fn = datadir + '/bbcar_%s_intindex.csv' % cf_name
        cf = pd.read_csv(cf_fn, index_col=0)
        cfte = cf.loc[test_indices,:]
        # using the best parameters, find the training record for best performing model 
        train_record_fn = train_record_dir + '/' + gen_feature + '_' + cf_name + '/scanmap_k%s_wcls%s_C%s_lr%s.p' % (nc, wcls, C, lr)
        try:
            record = torch.load(train_record_fn, map_location=torch.device('cpu'))
            Wte = record["state_dict"]['Wte']
            fc_weight = record["state_dict"]['fc.weight']
            fc_bias = record["state_dict"]['fc.bias']
            # get model output
            try:
                y_pred = np.round(softmax(np.dot(np.concatenate((cfte, Wte), axis=1), fc_weight.T) + np.array(fc_bias), axis=1)).astype(int)[:,1]
                np.savetxt(resultsdir + '/pred_%s_%s_k%s_wcls%s_C%s_lr%s.csv' % 
                    (gen_feature, cf_name, nc, wcls, C, lr), y_pred, delimiter=",")
            except ValueError:
                print('Genomic feature: %s' % gen_feature)
                print('Confounding features: %s' % cf_name)
                print('nc, C, wcls, lr = %s, %s, %s, %s' % (nc, C, wcls, lr))
            # y_pred.to_csv('bbcar/scanmap/results/pred_%s_%s_k%s_wcls%s_C%s_lr%s.csv' % 
            #     (gen_feature, cf_name, nc, wcls, C, lr))
        except FileNotFoundError:
            print('Could not load file: %s' % fn)
            continue
    except:
        print('Could not read file %s' % fn)

########################################################
#### Now save predictions of no-CNV baseline models ####
########################################################

from sklearn.linear_model import LogisticRegression

#### Load train, val, test indices #### 

train_ix = list(pd.read_csv(os.path.join(ixdir, 'train_indices_0.1val_0.2te.csv'), header=None, index_col=0).index)
val_ix = list(pd.read_csv(os.path.join(ixdir, 'val_indices_0.1val_0.2te.csv'), header=None, index_col=0).index)
test_ix = list(pd.read_csv(os.path.join(ixdir, 'test_indices_0.1val_0.2te.csv'), header=None, index_col=0).index)

y_train = pd.read_csv(os.path.join(datadir, 'bbcar_label_intindex.csv'), header=0, index_col=0).loc[train_ix,:].to_numpy().ravel()

results = pd.read_csv(os.path.join(resultsdir, 'baseline_additional_features.txt'))

for cf_name in results['features'].unique():
    subresults = results.iloc[results['features'].values==cf_name,:]
    bestpf = subresults.sort_values('val bal acc', ascending=False).iloc[0,:]
    C = bestpf['C']
    if C >=1:
        C = int(C)
    cf_fn = datadir + '/bbcar_%s_intindex.csv' % cf_name
    cf = pd.read_csv(cf_fn, index_col=0)
    cftr, cfte = cf.loc[train_ix,:], cf.loc[test_ix,:]
    lrm = LogisticRegression(C=C, max_iter=2000)
    lrm.fit(cftr, y_train)
    y_test_pred = np.array(lrm.predict(cfte)).reshape(-1,1)
    np.savetxt(resultsdir + '/pred_noCNV_%s_C%s.csv' % 
        (cf_name, C), y_test_pred, delimiter=",")
