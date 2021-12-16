import os
import sys
import numpy as np
import pandas as pd
import torch
import glob
import torch
import matplotlib.pyplot as plt

niter=2000
seed=1
outdir='/projects/b1042/ClareLab/saya/train_record_mlp'

for fn in glob.glob('bbcar/model_performance/results_mlp/results_*.txt'):
    try:
        # find best performing model's parameters
        results = pd.read_csv(fn, skiprows=1)
        bestpf = results.sort_values('val bal acc', ascending=False).iloc[0,:]
        C, lr = bestpf['C'], bestpf['lr']
        if C >=1:
            C = int(C)
        # get feature names 
        gen_feature = fn.split('.')[0].split('/')[-1].split('_')[1]
        # using the best parameters, find the training record for best performing model 
        train_record_fn = '%s/bbcarmlp_%s_C%s_lr%s_niter%s_seed%s.p' % (outdir, gen_feature, C, lr, niter, seed)
        try:
            record = torch.load(train_record_fn, map_location=torch.device('cpu'))
            plt.plot(record['report']['loss'])
            plt.xlabel('Epochs')
            plt.ylabel('Loss')
            plt.title('%s, C=%s, lr=%s' % (gen_feature, C, lr))
            plt.savefig('bbcar/model_performance/results_mlp/learning_curve_mlp_%s_C%s_lr%s.png' % (gen_feature, C, lr))
            plt.close()
        except FileNotFoundError:
            print('Could not load file: %s' % train_record_fn)
            continue
    except:
        print('Could not read file %s' % fn)
        continue
