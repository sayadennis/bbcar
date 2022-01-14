import os
import sys
import numpy as np
import pandas as pd
import torch
import glob
from scipy.special import softmax
import matplotlib.pyplot as plt

####################################
#### First plot learning curves ####
####################################

record_dn = '/projects/b1042/ClareLab/saya/train_record_siamese/regthres_clin'

## C = 0.001

for wcls in [0.01, 0.1, 0.5, 1, 2, 10, 50]:
    fn = os.path.join(record_dn, 'scanmap_k30_wcls%s_C0.001_lr1e-05.p' % wcls)
    record = torch.load(fn, map_location=torch.device('cpu'))
    plt.plot(record['report']['loss'], label='wcls=%s' % wcls)

plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()
plt.title('regthres, clin, k=30, C=0.001, lr=1e-05')
plt.savefig('bbcar/scanmap/results_siamese/learning_curve_scanmap_regthres_clin_k30_wclsrange_C0.001_lr1e-05.png')
plt.close()

## C = 0.01

for wcls in [0.01, 0.1, 0.5]:
    fn = os.path.join(record_dn, 'scanmap_k30_wcls%s_C0.01_lr1e-05.p' % wcls)
    record = torch.load(fn, map_location=torch.device('cpu'))
    plt.plot(record['report']['loss'], label='wcls=%s' % wcls)

plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()
plt.title('regthres, clin, k=30, C=0.001, lr=1e-05')
plt.savefig('bbcar/scanmap/results_siamese/learning_curve_scanmap_regthres_clin_k30_wclsrange_C0.01_lr1e-05.png')
plt.close()


##########################################
#### Next summarize model performance ####
##########################################

from sklearn.metrics import accuracy_score, balanced_accuracy_score, precision_score, recall_score, f1_score

scanmapdir = '/projects/b1122/saya/scanmap_data'

test_indices = list(pd.read_csv(os.path.join(scanmapdir, 'test_indices_0.1val_0.2te.csv'), header=None, index_col=0).index)
y_test = pd.read_csv(os.path.join(scanmapdir, 'bbcar_label_intindex.csv'), header=0, index_col=0).loc[test_indices,:].to_numpy().ravel()

pf_summary = pd.DataFrame(columns=['C', 'wcls', 'acc', 'bal acc', 'precision', 'recall', 'f1'])

nc = 30
lr = 1e-05
gen_feature = 'regthres'
cf_name = 'clin'
cf_fn = '/projects/b1122/saya/scanmap_data/bbcar_%s_intindex.csv' % cf_name
cf = pd.read_csv(cf_fn, index_col=0)
cfte = cf.loc[test_indices,:]

def siamese_predict(w, cf, fc_weight, fc_bias):
    # w = torch.from_numpy(w).float()
    cf = torch.from_numpy(np.array(cf)).float()
    # create index map and paired input matrices for W and cf 
    index_map = torch.cat([torch.repeat_interleave(torch.arange(w.shape[0]), w.shape[0], dim=0).reshape(1,-1), torch.arange(w.shape[0]).repeat(1, w.shape[0])], dim=0)
    w_pairs = torch.cat([torch.repeat_interleave(w, w.shape[0], dim=0), w.repeat(w.shape[0], 1)], dim=1)
    cf_pairs = torch.cat([torch.repeat_interleave(cf, cf.shape[0], dim=0), cf.repeat(cf.shape[0], 1)], dim=1)
    # create four-class prediction: list of values {0,1,2,3} 
    y_pairs_fourclass = torch.argmax(torch.mm(torch.cat([w_pairs, cf_pairs],1), fc_weight.T) + fc_bias, dim=1)
    # create n_samples x 2 target with values {0,1}
    y_pairs = torch.tensor([], dtype=int)
    for i in range(len(y_pairs_fourclass)):
        if torch.eq(y_pairs_fourclass[i],torch.tensor([0])):
            y_pairs = torch.cat([y_pairs, torch.tensor([0,0])], dim=0)
        elif torch.eq(y_pairs_fourclass[i],torch.tensor([1])):
            y_pairs = torch.cat([y_pairs, torch.tensor([1,0])], dim=0)
        elif torch.eq(y_pairs_fourclass[i],torch.tensor([2])):
            y_pairs = torch.cat([y_pairs, torch.tensor([0,1])], dim=0)
        elif torch.eq(y_pairs_fourclass[i],torch.tensor([3])):
            y_pairs = torch.cat([y_pairs, torch.tensor([1,1])], dim=0)
    y_pairs = y_pairs.reshape(-1,2)
    y_pairs = torch.t(y_pairs)
    # now make prediction vector y by taking "votes" from predictions on y_pairs
    y = torch.tensor([], dtype=int)
    for i in range(w.shape[0]):
        mapping = torch.where(index_map==i)
        votes = []
        for j,k in zip(mapping[0], mapping[1]):
            votes.append(y_pairs[j,k])
        if torch.sum(torch.tensor(votes))/len(votes) >= 0.5:
            y = torch.cat([y, torch.tensor([1])], dim=0)
        else:
            y = torch.cat([y, torch.tensor([0])], dim=0)
    return y


C = 0.001

for wcls in [0.01, 0.1, 0.5, 1, 2, 10, 50]:
    train_record_fn = '/projects/b1042/ClareLab/saya/train_record_siamese/' + gen_feature + '_' + cf_name + '/scanmap_k%s_wcls%s_C%s_lr%s.p' % (nc, wcls, C, lr)
    record = torch.load(train_record_fn, map_location=torch.device('cpu'))
    Wte = record["state_dict"]['Wte']
    fc_weight = record["state_dict"]['fc.weight']
    fc_bias = record["state_dict"]['fc.bias']
    # y_pred = np.round(softmax(np.dot(np.concatenate((cfte, Wte), axis=1), fc_weight.T) + np.array(fc_bias), axis=1)).astype(int)[:,1]
    y_pred = siamese_predict(Wte, cfte, fc_weight, fc_bias)
    pf_summary = pf_summary.append({
        'C' : C,
        'wcls' : wcls,
        'acc' : accuracy_score(y_test, y_pred),
        'bal acc' : balanced_accuracy_score(y_test, y_pred),
        'precision' : precision_score(y_test, y_pred),
        'recall' : recall_score(y_test, y_pred),
        'f1' : f1_score(y_test, y_pred)
    }, ignore_index=True)


C = 0.01

for wcls in [0.01, 0.1, 0.5]:
    train_record_fn = '/projects/b1042/ClareLab/saya/train_record_siamese/' + gen_feature + '_' + cf_name + '/scanmap_k%s_wcls%s_C%s_lr%s.p' % (nc, wcls, C, lr)
    record = torch.load(train_record_fn, map_location=torch.device('cpu'))
    Wte = record["state_dict"]['Wte']
    fc_weight = record["state_dict"]['fc.weight']
    fc_bias = record["state_dict"]['fc.bias']
    # y_pred = np.round(softmax(np.dot(np.concatenate((cfte, Wte), axis=1), fc_weight.T) + np.array(fc_bias), axis=1)).astype(int)[:,1]
    y_pred = siamese_predict(Wte, cfte, fc_weight, fc_bias)
    pf_summary = pf_summary.append({
        'C' : C,
        'wcls' : wcls,
        'acc' : accuracy_score(y_test, y_pred),
        'bal acc' : balanced_accuracy_score(y_test, y_pred),
        'precision' : precision_score(y_test, y_pred),
        'recall' : recall_score(y_test, y_pred),
        'f1' : f1_score(y_test, y_pred)
    }, ignore_index=True)
