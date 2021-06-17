import os
import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score

def get_pval(y_a, y_b, y_true, n_repeat=10000): # null hypothesis is that y_a is not better than y_b
    # make sure y_pred and y_true are numpy arrays 
    if type(y_a) != np.ndarray:
        y_a = y_a.to_numpy(dtype=int)
    else:
        y_a = np.array(y_a, dtype=int)
    if type(y_b) != np.ndarray:
        y_b = y_b.to_numpy(dtype=int)
    else:
        y_b = np.array(y_b, dtype=int)
    # convert to 1D array 
    y_a = y_a.ravel()
    y_b = y_b.ravel()

    pf_ya = balanced_accuracy_score(y_true, y_a)
    pf_yb = balanced_accuracy_score(y_true, y_b)
    
    # observed score 
    obs_diff = abs(pf_ya - pf_yb)
    
    # simulation 
    comb_array = np.array([y_a, y_b]).T
    diff_record = [] # performance difference record
    for i in range(n_repeat):
        sim_array = np.copy(comb_array)
        for j in range(comb_array.shape[0]):
            sim_array[j,:] = np.random.permutation(sim_array[j,:])
        # record performance
        diff_record.append(abs(balanced_accuracy_score(y_true, sim_array[:,0].ravel()) - balanced_accuracy_score(y_true, sim_array[:,1].ravel())))
    p = (np.sum(np.array(diff_record) >= obs_diff)+1)/(len(diff_record)+1)
    return p, pf_ya, pf_yb
