import os
import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, accuracy_score, precision_score, recall_score, f1_score

#### Load data and label #### 

dn = '/projects/b1122/saya/scanmap_data'

clin = pd.read_csv(os.path.join(dn, 'bbcar_clin_intindex.csv'), index_col=0)
mut = pd.read_csv(os.path.join(dn, 'bbcar_mutsig_intindex.csv'), index_col=0)
driv = pd.read_csv(os.path.join(dn, 'bbcar_driversomatic_intindex.csv'), index_col=0)
prs = pd.read_csv(os.path.join(dn, 'bbcar_prs_intindex.csv'), index_col=0)

clin_mut = pd.read_csv(os.path.join(dn, 'bbcar_clin_mut_intindex.csv'), index_col=0)
clin_driv = pd.read_csv(os.path.join(dn, 'bbcar_clin_driversomatic_intindex.csv'), index_col=0)
clin_prs = pd.read_csv(os.path.join(dn, 'bbcar_clin_prs_intindex.csv'), index_col=0)

clin_mut_prs = pd.read_csv(os.path.join(dn, 'bbcar_clin_mut_prs_intindex.csv'), index_col=0)
clin_mut_driv = pd.read_csv(os.path.join(dn, 'bbcar_clin_mut_driversomatic_intindex.csv'), index_col=0)
clin_prs_driv = pd.read_csv(os.path.join(dn, 'bbcar_clin_prs_driversomatic_intindex.csv'), index_col=0)

label = pd.read_csv(os.path.join(dn, 'bbcar_label_intindex.csv'), index_col=0)

#### Load train, val, test indices #### 

train_ix = list(pd.read_csv(os.path.join(dn, 'train_indices_0.1val_0.2te.csv'), header=None, index_col=0).index)
val_ix = list(pd.read_csv(os.path.join(dn, 'val_indices_0.1val_0.2te.csv'), header=None, index_col=0).index)
test_ix = list(pd.read_csv(os.path.join(dn, 'test_indices_0.1val_0.2te.csv'), header=None, index_col=0).index)

#### Train models and record performance #### 

comb_dict = {
    'clin' : clin,
    'mut' : mut,
    'driversomatic' : driv,
    'prs' : prs,
    'clin_mut' : clin_mut,
    'clin_driversomatic' : clin_driv,
    'clin_prs' : clin_prs,
    'clin_mut_prs' : clin_mut_prs,
    'clin_mut_driversomatic' : clin_mut_driv,
    'clin_prs_driversomatic' : clin_prs_driv
}

Cs = [0.001, 0.01, 0.1, 1, 10, 100, 1000]

print('features,C,tr bal acc,val bal acc,te acc,te bal acc,te precis,te recall,te f1')
for key in comb_dict.keys():
    X_train, X_val, X_test = comb_dict[key].iloc[train_ix], comb_dict[key].iloc[val_ix], comb_dict[key].iloc[test_ix]
    y_train, y_val, y_test = label.iloc[train_ix].values.ravel(), label.iloc[val_ix].values.ravel(), label.iloc[test_ix].values.ravel()
    for C in Cs:
        lrm = LogisticRegression(C=C, max_iter=2000)
        lrm.fit(X_train, y_train)
        # record training performance
        y_train_pred = lrm.predict(X_train)
        acctr = balanced_accuracy_score(y_train, y_train_pred)
        # record validation performance
        y_val_pred = lrm.predict(X_val)
        accval = balanced_accuracy_score(y_val, y_val_pred)
        # record test set performance
        y_test_pred = np.array(lrm.predict(X_test)).reshape(-1,1)
        np.savetxt('bbcar/scanmap/results/pred_noCNV_%s_C%s.csv' % 
            (key, C), y_test_pred, delimiter=",")
        accte = accuracy_score(y_test, y_test_pred)
        balaccte = balanced_accuracy_score(y_test, y_test_pred)
        preciste = precision_score(y_test, y_test_pred)
        recallte = recall_score(y_test, y_test_pred)
        f1te = f1_score(y_test, y_test_pred)
        print('%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f' % 
            (key, C, acctr, accval, accte, balaccte, preciste, recallte, f1te)) # best_iter, mse, mse_tr, mse_val, mse_te
