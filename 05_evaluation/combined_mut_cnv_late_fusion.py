import pickle
from collections import defaultdict 

import numpy as np
import pandas as pd
from sklearn import metrics
from sklearn.model_selection import StratifiedKFold

project_dn = '/projects/b1131/saya/bbcar'

###############################################
#### Models from original feature matrices ####
###############################################

mut_dn = f'{project_dn}/data/02a_mutation/08_feature_matrix/20230423_signature_results'
mut_fn = 'sbs_96_original_per_sample.csv'
cnv_dn = f'{project_dn}/data/02b_cnv/inhouse_signatures'
cnv_fn = f'inhouse_sig_batcheffect_rm_combat.csv'
target_dn = '/projects/b1131/saya/bbcar/data/clinical'
target_fn = 'bbcar_label_studyid_from_gatk_filenames.csv'
indexdir = '/projects/b1131/saya/bbcar/train_test_splits'

mut_modeldir = f'{project_dn}/model_interpretations/202310_original_96'
cnv_modeldir = f'{project_dn}/model_interpretations/cnv_sig_nonNMF/inhouse_sig_batcheffect_rm_combat'

mut_model_fn = '20231010_saved_best_XGB_sbs_96_original_per_sample.p'
cnv_model_fn = '20231011_saved_best_XGB_inhouse_sig_batcheffect_rm_combat.p'

with open(f'{mut_modeldir}/{mut_model_fn}', 'rb') as f:
    mut_model = pickle.load(f)

with open(f'{cnv_modeldir}/{cnv_model_fn}', 'rb') as f:
    cnv_model = pickle.load(f)

X_mut = pd.read_csv(f'{mut_dn}/{mut_fn}', index_col=0)
X_cnv = pd.read_csv(f'{cnv_dn}/{cnv_fn}', index_col=0)
y = pd.read_csv(f'{target_dn}/{target_fn}', index_col=0)

train_ix = pd.read_csv(f'{indexdir}/train_index.txt', header=None).to_numpy().ravel()
test_ix = pd.read_csv(f'{indexdir}/test_index.txt', header=None).to_numpy().ravel()

# Split data
X_mut_train, X_mut_test = X_mut.loc[train_ix,:], X_mut.loc[test_ix,:]
X_cnv_train, X_cnv_test = X_cnv.loc[train_ix,:], X_cnv.loc[test_ix,:]
y_train, y_test = y.loc[train_ix,:], y.loc[test_ix,:]

## Simple average of the two scores 
def late_fusion_model(X_mut, X_cnv, mut_model, cnv_model, weight=0.5):
    mut_proba = mut_model.predict_proba(X_mut)[:,1]
    cnv_proba = cnv_model.predict_proba(X_cnv)[:,1]
    prob = np.average([mut_proba, cnv_proba], axis=0, weights=[(1-weight), weight])
    pred = np.array(prob>0.5, dtype=float)
    return pred, prob

y_train_pred, y_train_prob = late_fusion_model(X_mut_train, X_cnv_train, mut_model, cnv_model)
y_test_pred, y_test_prob = late_fusion_model(X_mut_test, X_cnv_test, mut_model, cnv_model)

def evaluate_model(
        y_true: np.ndarray, y_pred: np.ndarray, y_prob: np.ndarray, 
    ):
    performance = pd.DataFrame(
        index=['Performance'],
        columns=['Bal Acc', 'ROC AUC', 'Precision', 'Recall', 'F1'],    
    )   
    performance.loc['Performance', 'Bal Acc'] = metrics.balanced_accuracy_score(y_true, y_pred)
    performance.loc['Performance', 'ROC AUC'] = metrics.roc_auc_score(y_true, y_prob)
    performance.loc['Performance', 'Precision'] = metrics.precision_score(y_true, y_pred)
    performance.loc['Performance', 'Recall'] = metrics.recall_score(y_true, y_pred)
    performance.loc['Performance', 'F1'] = metrics.f1_score(y_true, y_pred)
    print(performance.T)

print(f'\n\n#### Original feature matrix & Simple average ####')
print(f'## Train set performance ##')
evaluate_model(y_train, y_train_pred, y_train_prob)

print(f'## Test set performance ##')
evaluate_model(y_test, y_test_pred, y_test_prob)

## Weighted average of the two scores
skf = StratifiedKFold(n_splits=5, random_state=42, shuffle=True)
weights = np.arange(0.1, 1.0, step=0.1)
scores = defaultdict(list)
for weight in weights:
    for i, (train_index, test_index) in enumerate(skf.split(X_mut_train, y_train)):
        # Split dataset
        X_mut_train_folds = X_mut_train.iloc[train_index,:]
        X_cnv_train_folds = X_cnv_train.iloc[train_index,:]
        y_train_folds = y_train.iloc[train_index,:]
        X_mut_test_fold = X_mut_train.iloc[test_index,:]
        X_cnv_test_fold = X_cnv_train.iloc[test_index,:]
        y_test_fold = y_train.iloc[test_index,:]
        # Calculate and record score for this fold
        y_pred_fold, y_score_fold = late_fusion_model(X_mut_test_fold, X_cnv_test_fold, mut_model, cnv_model, weight=weight)
        i_fold_score = metrics.roc_auc_score(y_test_fold, y_score_fold)
        scores[weight].append(i_fold_score)

print(f'\n\n#### Original feature matrix & Weighted average ####')

for weight in weights:
    print(f'Mean score across folds for weight {weight:.1f}: {np.average(scores[weight])}')
    y_test_pred, y_test_score = late_fusion_model(X_mut_test, X_cnv_test, mut_model, cnv_model, weight=weight)
    test_score = metrics.roc_auc_score(y_test, y_test_score)
    print(f'Test set score with weight {weight:.1f}: {test_score}\n')

########################################
#### Models from signature matrices ####
########################################

mut_dir = '202310_signature'
cnv_dir = 'cnv_sig_nonNMF/inhouse_cnv_sig_per_sample'

## Simple average of the two scores 

## Weighted average of the two scores 


