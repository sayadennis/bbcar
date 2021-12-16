import os
import sys
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn import metrics

inputfn = '/projects/b1122/saya/06_modified_data/cyto_copy_conf90_studyindex.csv'
labfn = '/projects/b1122/saya/bbcar_label_studyindex.csv'
gisticdir = '/projects/b1122/saya/03_gistic2_out_conf90/all'
indexdir='/projects/b1122/saya/indices'

X = pd.read_csv(inputfn, index_col=0)
y = pd.read_csv(labfn, index_col=0)

train_ix = pd.read_csv('%s/train_ix.csv' % (indexdir), header=None)
test_ix = pd.read_csv('%s/test_ix.csv' % (indexdir), header=None)

X_train, X_test = X.iloc[train_ix[0]], X.iloc[test_ix[0]]
y_train, y_test = y.iloc[train_ix[0]], y.iloc[test_ix[0]]

X_train = X_train.to_numpy(dtype=int)
X_test = X_test.to_numpy(dtype=int)

y_train = y_train.to_numpy(dtype=int)
y_train = y_train.ravel()
y_test = y_test.to_numpy(dtype=int)
y_test = y_test.ravel()

lasso = LogisticRegression(C=10., penalty="l1", solver="liblinear", class_weight="balanced", max_iter=1000, random_state=0)

lasso.fit(X_train, y_train)

y_pred = lasso.predict(X_test)
y_prob = lasso.predict_proba(X_test)[:,1]

coefs = pd.DataFrame(lasso.coef_[0][np.where(lasso.coef_[0]!=0)], index=X.columns[np.where(lasso.coef_[0]!=0)], columns=['coefs'])

coefs.sort_values('coefs', ascending=False).to_csv('bbcar/cyto_cnv_lasso_nonzero_coefs.csv', header=True, index=True)

################################################################################
#### Get a list of genes that are in the amp/del regions of these cytobands ####
################################################################################

coefs = pd.read_csv('bbcar/cyto_cnv_lasso_nonzero_coefs.csv', index_col=0)

amp_genes = pd.read_csv(f'{gisticdir}/amp_genes.conf_90.txt', sep='\t', index_col=0)
del_genes = pd.read_csv(f'{gisticdir}/del_genes.conf_90.txt', sep='\t', index_col=0)

amp_genes = amp_genes.iloc[3:,[x in coefs.index for x in amp_genes.columns]].values.ravel()
del_genes = del_genes.iloc[3:,[x in coefs.index for x in del_genes.columns]].values.ravel()

amp_genes = list(amp_genes[~pd.isnull(amp_genes)])
del_genes = list(del_genes[~pd.isnull(del_genes)])

with open('bbcar/lasso_nonzero_amp_genes.txt', 'w') as f:
    for item in amp_genes:
        f.write(item + '\n')

with open('bbcar/lasso_nonzero_del_genes.txt', 'w') as f:
    for item in del_genes:
        f.write(item + '\n')

