import pickle
import numpy as np
import pandas as pd
from sklearn.metrics import confusion_matrix, roc_curve, auc, roc_auc_score
import matplotlib.pyplot as plt

proj_dn = '/projects/b1131/saya/bbcar'

inputdir = f'{proj_dn}/data/02a_mutation/08_feature_matrix/'
labeldir = f'{proj_dn}/data/clinical/'
outdir = f'{proj_dn}/model_interpretations'
ixdir = f'{proj_dn}/train_test_splits'

###################
#### Prep data ####
###################

# Load performance 
pf = pd.read_csv(f'{outdir}/bbcar_genesommut_roc_auc.csv', index_col=0)

# Load model
with open(f'{proj_dn}/models/20221102_saved_best_SVM_cts_per_gene_Polyphen2_D.p', 'rb') as f:
    m = pickle.load(f)

# Load data 
X = pd.read_csv(f'{inputdir}/cts_per_gene_Polyphen2_D.csv', index_col=0)
y = pd.read_csv(f'{labeldir}/bbcar_redcap_label_studyid.csv', index_col=0)

# Split data 
train_ix = pd.read_csv(f'{ixdir}/train_ix.csv', header=None).to_numpy().ravel()
test_ix = pd.read_csv(f'{ixdir}/test_ix.csv', header=None).to_numpy().ravel()

X_train, X_test = X.iloc[[i in train_ix for i in X.index],:], X.iloc[[i in test_ix for i in X.index],:]
y_train, y_test = y.iloc[[i in train_ix for i in y.index],:], y.iloc[[i in test_ix for i in y.index],:]

# Align order
y_train = y_train.loc[X_train.index,:]
y_test = y_test.loc[X_test.index,:]

####################################
#### Get predictions and scores ####
####################################

y_test_pred = m.predict(X_test)
y_test_score = m.predict_proba(X_test)[:,1]

##############################
#### Get confusion matrix ####
##############################

pd.DataFrame(
    confusion_matrix(y_test, y_test_pred),
    index=['true neg', 'true pos'],
    columns=['pred neg', 'pred pos']
).to_csv(f'{proj_dn}/out/confus_mx_gene_wise_mut_cts.csv')

########################
#### Plot ROC curve ####
########################

fpr, tpr, _ = roc_curve(y_test, y_test_score)
roc_auc = auc(fpr, tpr)

plt.figure()
lw = 2
plt.plot(
    fpr,
    tpr,
    color='darkorange',
    lw=lw,
    label='ROC curve (area = %0.2f)' % roc_auc,
)
plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('ROC for gene-wise mutation counts')
plt.legend(loc='lower right')
plt.savefig(f'{proj_dn}/plots/gene_wise_mut_cts_roc.png')
plt.close()
