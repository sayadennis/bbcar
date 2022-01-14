import os
import numpy as np
import pandas as pd
import glob

def mx_sumstats(mx):
    """
    This function takes a feature matrix (Pandas DataFrame) and calculates, across all entries, their summary statistics.
    It returns a table with all the summary stats 
    """
    sumstat_dict = {}
    sumfeatures = ["# features", "Mean", "Min", "25%", "Median", "75%", "Max", "% positive entries", "% negative entries", "% zero entries"]
    for sumfeature in sumfeatures:
        sumstat_dict[sumfeatures[0]] = len(mx.columns) # number of features 
        sumstat_dict[sumfeatures[1]] = np.mean(mx.to_numpy().ravel()) # mean of all entries 
        sumstat_dict[sumfeatures[2]] = np.quantile(mx.to_numpy().ravel(), 0)
        sumstat_dict[sumfeatures[3]] = np.quantile(mx.to_numpy().ravel(), 0.25)
        sumstat_dict[sumfeatures[4]] = np.quantile(mx.to_numpy().ravel(), 0.5)
        sumstat_dict[sumfeatures[5]] = np.quantile(mx.to_numpy().ravel(), 0.75)
        sumstat_dict[sumfeatures[6]] = np.quantile(mx.to_numpy().ravel(), 1)
        sumstat_dict[sumfeatures[7]] = 100*np.sum(mx.to_numpy().ravel()>0)/len(mx.to_numpy().ravel())
        sumstat_dict[sumfeatures[8]] = 100*np.sum(mx.to_numpy().ravel()<0)/len(mx.to_numpy().ravel())
        sumstat_dict[sumfeatures[9]] = 100*np.sum(mx.to_numpy().ravel()==0)/len(mx.to_numpy().ravel())
    return sumstat_dict

def applyThres(table, thres=0.05):
    # table: Pandas DataFrame
    # thres: float
    bool_col = [] # list of boolean indicating the columns to keep at True
    for colname in table.columns:
        ct = 0
        for i in range(len(table[colname])):
            if table[colname][i] != 0:
                ct += 1
            else:
                continue
        rate = ct/len(table[colname])
        if ((rate < thres) | (rate > (1-thres))): # change this so that we only eliminate features that are rare?? (because common can add value?)
            bool_col.append(False)
        else:
            bool_col.append(True)    
    return table.iloc[:,bool_col]

def applyThres_featuresize(table, thres=[0.05, 0.10, 0.20, 0.25]):
    # fs_table is a table summarizing how feature size is affected by freq-based reduction 
    fs_table = pd.DataFrame(None, columns=["Feature size"], index=["0.05", "0.10", "0.20", "0.25"])
    fs_table.loc["Original"]["Feature size"] = table.shape[1]
    fs_table.loc["0.05"]["Feature size"] = applyThres(table, thres=0.05).shape[1]
    fs_table.loc["0.10"]["Feature size"] = applyThres(table, thres=0.10).shape[1]
    fs_table.loc["0.20"]["Feature size"] = applyThres(table, thres=0.20).shape[1]
    fs_table.loc["0.25"]["Feature size"] = applyThres(table, thres=0.25).shape[1]
    return fs_table
