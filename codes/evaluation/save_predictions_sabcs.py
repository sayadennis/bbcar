import os
import sys
import numpy as np
import pandas as pd
from scipy.special import softmax
import glob
import torch

def create_pairs(x): # x is a numpy array 
    index_map = torch.cat([
        torch.repeat_interleave(torch.arange(x.shape[0]), x.shape[0], dim=0).reshape(1,-1), 
        torch.arange(x.shape[0]).repeat(1, x.shape[0])
    ], dim=0)
    if type(x)!=torch.Tensor:
        x = torch.from_numpy(np.array(x)).float()
    x_pairs = torch.cat([torch.repeat_interleave(x, x.shape[0], dim=0), x.repeat(x.shape[0], 1)], dim=1)
    return index_map, np.array(x_pairs)

def predict_siamese(fc_weight, fc_bias, Wte, cfte=None):
    index_map, Wte_pairs = create_pairs(Wte)
    if cfte is not None:
        _, cfte_pairs = create_pairs(cfte)
    # get one-hot vector 
    if cfte is not None:
        y_pred_fourclass = torch.from_numpy(np.argmax(softmax(np.dot(np.concatenate((Wte_pairs, cfte_pairs), axis=1), fc_weight.T) + np.array(fc_bias), axis=1), axis=1))
    else:
        y_pred_fourclass = torch.from_numpy(np.argmax(softmax(np.dot(Wte_pairs, fc_weight.T) + np.array(fc_bias), axis=1), axis=1))
    # create the pairs predictions based on the one-hot vector 
    y_pairs = torch.tensor([], dtype=int)
    for i in range(y_pred_fourclass.shape[0]):
        if torch.eq(y_pred_fourclass[i],torch.tensor([0])):
            y_pairs = torch.cat([y_pairs, torch.tensor([0,0])], dim=0)
        elif torch.eq(y_pred_fourclass[i],torch.tensor([1])):
            y_pairs = torch.cat([y_pairs, torch.tensor([1,0])], dim=0)
        elif torch.eq(y_pred_fourclass[i],torch.tensor([2])):
            y_pairs = torch.cat([y_pairs, torch.tensor([0,1])], dim=0)
        elif torch.eq(y_pred_fourclass[i],torch.tensor([3])):
            y_pairs = torch.cat([y_pairs, torch.tensor([1,1])], dim=0)
    y_pairs = y_pairs.reshape(-1,2)
    y_pairs = torch.t(y_pairs)
    # take the vote to get prediction for each sample 
    y = torch.tensor([], dtype=int)
    for i in range(Wte.shape[0]):
        mapping = torch.where(index_map==i)
        votes = []
        for j,k in zip(mapping[0], mapping[1]):
            votes.append(y_pairs[j,k])
        if torch.sum(torch.tensor(votes))/len(votes) >= 0.5:
            y = torch.cat([y, torch.tensor([1])], dim=0)
        else:
            y = torch.cat([y, torch.tensor([0])], dim=0)
    # np.savetxt(outfn, y, delimiter=",")
    return y


datadir = '/projects/b1122/saya/06_modified_data'
cfdir = '/projects/b1122/saya/bbcar_non_cnv_features'
labdir = '/projects/b1122/saya'
ixdir = '/projects/b1122/saya/indices'
resultsdir = 'bbcar/model_performance'
train_record_dir = '/projects/b1042/ClareLab/saya/train_record_siamese'

test_indices = list(pd.read_csv(os.path.join(ixdir, 'test_indices_0.1val_0.2te.csv'), header=None, index_col=0).index)
y_test = pd.read_csv(os.path.join(labdir, 'bbcar_label_intindex.csv'), header=0, index_col=0).loc[test_indices,:].to_numpy().ravel()

##################
#### CNA only ####
##################

gen_feature = 'regthres'
# cf_name = 'clin'
# cf_fn = cfdir + '/bbcar_%s_intindex.csv' % cf_name
# cf = pd.read_csv(cf_fn, index_col=0)
# cfte = cf.loc[test_indices,:]

fn = f'{resultsdir}/results_scanmap_siamese_nocf/results_scanmap_siamese_regthres_nocf.txt'
results = pd.read_csv(fn, skiprows=2)
bestpf = results.sort_values('val bal acc', ascending=False).iloc[0,:]
nc, C, wcls, lr = int(bestpf['nc']), bestpf['C'], bestpf['wcls'], bestpf['lr']

if wcls >= 1: 
    wcls = int(wcls)

if C >=1:
    C = int(C)

# using the best parameters, find the training record for best performing model 
train_record_fn = train_record_dir + '_nocf/scanmap_k%s_wcls%s_C%s_lr%s.p' % (nc, wcls, C, lr)
record = torch.load(train_record_fn, map_location=torch.device('cpu'))

# get factorized test matrix as well as weights and biases for the fully connected layer 
Wte = record["state_dict"]['Wte']
fc_weight = record["state_dict"]['fc.weight']
fc_bias = record["state_dict"]['fc.bias']

# get predicted results and save
y_pred = predict_siamese(fc_weight, fc_bias, Wte)
np.savetxt(f'{resultsdir}/results_scanmap_siamese_nocf/pred_{gen_feature}_nocf_k{nc}_wcls{wcls}_C{C}_lr{lr}.txt', y_pred.numpy())

########################
#### CNA + clinical ####
########################

gen_feature = 'regthres'
cf_name = 'clin'
cf_fn = cfdir + '/bbcar_%s_intindex.csv' % cf_name
cf = pd.read_csv(cf_fn, index_col=0)
cfte = cf.loc[test_indices,:]

fn = f'{resultsdir}/results_scanmap_siamese/results_siamese_regthres_clin.txt'
results = pd.read_csv(fn, skiprows=2)
bestpf = results.sort_values('val bal acc', ascending=False).iloc[0,:]
nc, C, wcls, lr = int(bestpf['nc']), bestpf['C'], bestpf['wcls'], bestpf['lr']

if wcls >= 1: 
    wcls = int(wcls)

if C >=1:
    C = int(C)

# using the best parameters, find the training record for best performing model 
train_record_fn = train_record_dir + '/' + gen_feature + '_' + cf_name + '/scanmap_k%s_wcls%s_C%s_lr%s.p' % (nc, wcls, C, lr)
record = torch.load(train_record_fn, map_location=torch.device('cpu'))

# get factorized test matrix as well as weights and biases for the fully connected layer 
Wte = record["state_dict"]['Wte']
fc_weight = record["state_dict"]['fc.weight']
fc_bias = record["state_dict"]['fc.bias']

# get predicted results and save
y_pred = predict_siamese(fc_weight, fc_bias, Wte, cfte)
np.savetxt(f'{resultsdir}/results_scanmap_siamese/pred_{gen_feature}_{cf_name}_k{nc}_wcls{wcls}_C{C}_lr{lr}.txt', y_pred.numpy())

####################################
#### CNA + mutational signature ####
####################################

gen_feature = 'regthres'
cf_name = 'mut'
cf_fn = cfdir + '/bbcar_%s_intindex.csv' % cf_name
cf = pd.read_csv(cf_fn, index_col=0)
cfte = cf.loc[test_indices,:]

fn = f'{resultsdir}/results_scanmap_siamese/results_siamese_regthres_mut.txt'
results = pd.read_csv(fn, skiprows=2)
bestpf = results.sort_values('val bal acc', ascending=False).iloc[0,:]
nc, C, wcls, lr = int(bestpf['nc']), bestpf['C'], bestpf['wcls'], bestpf['lr']

if wcls >= 1: 
    wcls = int(wcls)

if C >=1:
    C = int(C)

# using the best parameters, find the training record for best performing model 
train_record_fn = train_record_dir + '/' + gen_feature + '_' + cf_name + '/scanmap_k%s_wcls%s_C%s_lr%s.p' % (nc, wcls, C, lr)
record = torch.load(train_record_fn, map_location=torch.device('cpu'))

# get factorized test matrix as well as weights and biases for the fully connected layer 
Wte = record["state_dict"]['Wte']
fc_weight = record["state_dict"]['fc.weight']
fc_bias = record["state_dict"]['fc.bias']

# get predicted results and save
y_pred = predict_siamese(fc_weight, fc_bias, Wte, cfte)
np.savetxt(f'{resultsdir}/results_scanmap_siamese/pred_{gen_feature}_{cf_name}_k{nc}_wcls{wcls}_C{C}_lr{lr}.txt', y_pred.numpy())

###############################################
#### CNA + clinical + mutational signature ####
###############################################

gen_feature = 'regthres'
cf_name = 'clin_mut'
cf_fn = cfdir + '/bbcar_%s_intindex.csv' % cf_name
cf = pd.read_csv(cf_fn, index_col=0)
cfte = cf.loc[test_indices,:]

fn = f'{resultsdir}/results_scanmap_siamese/results_regthres_clin_mut.txt'
results = pd.read_csv(fn, skiprows=2)
bestpf = results.sort_values('val bal acc', ascending=False).iloc[0,:]
nc, C, wcls, lr = int(bestpf['nc']), bestpf['C'], bestpf['wcls'], bestpf['lr']

if wcls >= 1: 
    wcls = int(wcls)

if C >=1:
    C = int(C)

# using the best parameters, find the training record for best performing model 
train_record_fn = train_record_dir + '/' + gen_feature + '_' + cf_name + '/scanmap_k%s_wcls%s_C%s_lr%s.p' % (nc, wcls, C, lr)
record = torch.load(train_record_fn, map_location=torch.device('cpu'))

# get factorized test matrix as well as weights and biases for the fully connected layer 
Wte = record["state_dict"]['Wte']
fc_weight = record["state_dict"]['fc.weight']
fc_bias = record["state_dict"]['fc.bias']

# get predicted results and save
y_pred = predict_siamese(fc_weight, fc_bias, Wte, cfte)
np.savetxt(f'{resultsdir}/results_scanmap_siamese/pred_{gen_feature}_{cf_name}_k{nc}_wcls{wcls}_C{C}_lr{lr}.txt', y_pred.numpy())

#####################################################
#### CNA + clinical + mutational signature + PRS ####
#####################################################

gen_feature = 'regthres'
cf_name = 'clin_mut_prs'
cf_fn = cfdir + '/bbcar_%s_intindex.csv' % cf_name
cf = pd.read_csv(cf_fn, index_col=0)
cfte = cf.loc[test_indices,:]

fn = f'{resultsdir}/results_scanmap_siamese/results_regthres_clin_mut_prs.txt'
results = pd.read_csv(fn, skiprows=2)
bestpf = results.sort_values('val bal acc', ascending=False).iloc[0,:]
nc, C, wcls, lr = int(bestpf['nc']), bestpf['C'], bestpf['wcls'], bestpf['lr']

if wcls >= 1: 
    wcls = int(wcls)

if C >=1:
    C = int(C)

# using the best parameters, find the training record for best performing model 
train_record_fn = train_record_dir + '/' + gen_feature + '_' + cf_name + '/scanmap_k%s_wcls%s_C%s_lr%s.p' % (nc, wcls, C, lr)
record = torch.load(train_record_fn, map_location=torch.device('cpu'))

# get factorized test matrix as well as weights and biases for the fully connected layer 
Wte = record["state_dict"]['Wte']
fc_weight = record["state_dict"]['fc.weight']
fc_bias = record["state_dict"]['fc.bias']

# get predicted results and save
y_pred = predict_siamese(fc_weight, fc_bias, Wte, cfte)
np.savetxt(f'{resultsdir}/results_scanmap_siamese/pred_{gen_feature}_{cf_name}_k{nc}_wcls{wcls}_C{C}_lr{lr}.txt', y_pred.numpy())

#########################################################################
#### CNA + clinical + mutational signature + driver somatic mutation ####
#########################################################################

gen_feature = 'regthres'
cf_name = 'clin_mut_driversomatic'
cf_fn = cfdir + '/bbcar_%s_intindex.csv' % cf_name
cf = pd.read_csv(cf_fn, index_col=0)
cfte = cf.loc[test_indices,:]

fn = f'{resultsdir}/results_scanmap_siamese/results_regthres_clin_mut_driversomatic.txt'
results = pd.read_csv(fn, skiprows=2)
bestpf = results.sort_values('val bal acc', ascending=False).iloc[0,:]
nc, C, wcls, lr = int(bestpf['nc']), bestpf['C'], bestpf['wcls'], bestpf['lr']

if wcls >= 1: 
    wcls = int(wcls)

if C >=1:
    C = int(C)

# using the best parameters, find the training record for best performing model 
train_record_fn = train_record_dir + '/' + gen_feature + '_' + cf_name + '/scanmap_k%s_wcls%s_C%s_lr%s.p' % (nc, wcls, C, lr)
record = torch.load(train_record_fn, map_location=torch.device('cpu'))

# get factorized test matrix as well as weights and biases for the fully connected layer 
Wte = record["state_dict"]['Wte']
fc_weight = record["state_dict"]['fc.weight']
fc_bias = record["state_dict"]['fc.bias']

# get predicted results and save
y_pred = predict_siamese(fc_weight, fc_bias, Wte, cfte)
np.savetxt(f'{resultsdir}/results_scanmap_siamese/pred_{gen_feature}_{cf_name}_k{nc}_wcls{wcls}_C{C}_lr{lr}.txt', y_pred.numpy())

