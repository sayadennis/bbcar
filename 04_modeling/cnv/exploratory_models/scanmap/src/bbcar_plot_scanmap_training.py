import os
import sys
import numpy as np
import pandas as pd
import torch
import glob
import torch
import matplotlib.pyplot as plt

for fn in glob.glob('bbcar/scanmap/results/results_*.txt'):
    try:
        # find best performing model's parameters
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
        # cf_fn = '/projects/b1122/saya/scanmap_data/bbcar_%s_intindex.csv' % cf_name
        # cf = pd.read_csv(cf_fn, index_col=0)
        # using the best parameters, find the training record for best performing model 
        train_record_fn = '/projects/b1042/ClareLab/saya/train_record/' + gen_feature + '_' + cf_name + '/scanmap_k%s_wcls%s_C%s_lr%s.p' % (nc, wcls, C, lr)
        try:
            record = torch.load(train_record_fn, map_location=torch.device('cpu'))
            plt.plot(record['report']['loss'])
            plt.xlabel('Epochs')
            plt.ylabel('Loss')
            plt.title('%s, %s, k=%d, wcls=%s, C=%s, lr=%s' % (gen_feature, cf_name, nc, wcls, C, lr))
            plt.savefig('bbcar/scanmap/results/learning_curve_scanmap_%s_%s_k%d_wcls%s_C%s_lr%s.png' % (gen_feature, cf_name, nc, wcls, C, lr))
            plt.close()
        except FileNotFoundError:
            print('Could not load file: %s' % train_record_fn)
            continue
    except:
        print('Could not read file %s' % fn)
        continue
