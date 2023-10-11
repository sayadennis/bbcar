import os
import sys
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt

from sklearn.metrics import accuracy_score, balanced_accuracy_score, roc_auc_score, precision_score, recall_score, f1_score, confusion_matrix, roc_curve

sys.path.append('classical_ml')
from ClassicalML import ClassicalML

din='/projects/b1122/saya/04_cleaned_cnv/all' # dir where input is stored
dlab='/projects/b1122/saya' # dir where label is stored
dix='/projects/b1122/saya/indices' # dir where indiced are stored
dpf='bbcar/model_performance/results_classicalml/all' # dir where model performance is stored
dout='bbcar/model_performance/results_classicalml/all' # dir where outputs of this script are saved

m_fn = '20210910_saved_best_LASSO_cyto_copy_conf90_all.p' # model iflename
pf_fn = 'bbcar_cytocopy_roc_auc.csv'

X = pd.read_csv(f'{din}/cyto_copy_conf90_all.csv', index_col=0)
y = pd.read_csv(f'{dlab}/bbcar_label_studyid.csv', index_col=0)

train_ix = pd.read_csv(f'{dix}/train_ix.csv', header=None)
test_ix = pd.read_csv(f'{dix}/test_ix.csv', header=None)
# train_ix = pd.read_csv(f'{dix}/train_studyid.csv', header=None) ## TRY USING THIS TOO ? 
# test_ix = pd.read_csv(f'{dix}/test_studyid.csv', header=None)

pf = pd.read_csv(f'{dpf}/{pf_fn}', index_col=0)

X_train, X_test = X.iloc[train_ix[0]], X.iloc[test_ix[0]]
y_train, y_test = y.iloc[train_ix[0]], y.iloc[test_ix[0]]

with open(m_fn, 'rb') as f:
    m = pickle.load(f)

# m.cv_results_
# m.best_estimator_
# m.best_score_
# m.scoring

## make sure that the metrics are correct 
y_pred = m.best_estimator_.predict(X_test)
y_prob = m.best_estimator_.predict_proba(X_test)[:,1]

confusion_matrix(y_test, y_pred, normalize=True)

accuracy_score(y_test,y_pred)
balanced_accuracy_score(y_test,y_pred)
roc_auc_score(y_test,y_prob)
precision_score(y_test,y_pred)
recall_score(y_test,y_pred)
f1_score(y_test,y_pred)

## inspect coefficients
lasso_coefs = pd.DataFrame({
    'cytoband' : X.columns[np.where(m.best_estimator_.coef_!=0)[1]],
    'lasso_coefs' : m.best_estimator_.coef_.ravel()[np.where(m.best_estimator_.coef_!=0)[1]]
}).sort_values('lasso_coefs', ascending=False)

coef_fn = m_fn.split('.')[0] + '_coefs.csv'
lasso_coefs.to_csv(f'{dout}/{coef_fn}', index=False, header=True)

## plot ROC 
plt.figure()
lw = 2
plt.plot(
    roc_curve(y_test, y_prob)[0],
    roc_curve(y_test, y_prob)[1],
    color="darkorange",
    lw=lw,
    label="ROC curve (area = %0.2f)" % roc_auc_score(y_test, y_prob),
)
plt.plot([0, 1], [0, 1], color="navy", lw=lw, linestyle="--")
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Receiver operating characteristic example")
plt.legend(loc="lower right")
plt.savefig('ROC_20210910_saved_best_LASSO_cyto_copy_conf90_all.png')
plt.close()
