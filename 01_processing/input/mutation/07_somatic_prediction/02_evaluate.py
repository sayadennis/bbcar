import os
import sys
import numpy as np
import pandas as pd
import pickle

# tools for evaluation
from sklearn import metrics
from numpy import interp

# visualization
import matplotlib.pyplot as plt

####################
#### Load stuff ####
####################

pon_source = sys.argv[1]

## Data
din = '/projects/b1131/saya/bbcar/data/02a_mutation/04_ml_features'
dout = f'/projects/b1131/saya/bbcar/data/02a_mutation/06_ml_eval/{pon_source}'

X = pd.read_csv(f'{din}/input_matched_{pon_source}.csv', index_col=0)
y = pd.read_csv(f'{din}/target_matched_{pon_source}.csv', index_col=0)

meta = pd.read_csv(f'{din}/04_imputed/features_imputed_{pon_source}.csv') # this contains meta information about each variant (variant exonic function etc.) 

## Model 
mfn = f'/projects/b1131/saya/bbcar/model_interpretations/{pon_source}/20230419_saved_best_XGB_input_matched_bbcar.p'
with open(mfn, 'rb') as f:
    m = pickle.load(f)

## Indices
train_ix = pd.read_csv(f'{din}/somatic_pred_ix/{pon_source}/train_index.txt', header=None)
test_ix = pd.read_csv(f'{din}/somatic_pred_ix/{pon_source}/test_index.txt', header=None)

X_train, X_test = X.loc[train_ix[0]], X.loc[test_ix[0]]
y_train, y_test = y.loc[train_ix[0]], y.loc[test_ix[0]]

##################
#### Evaluate ####
##################

y_train_pred = m.predict(X_train.to_numpy())
y_test_pred = m.predict(X_test.to_numpy())
y_test_prob = m.predict_proba(X_test.to_numpy())[:,1]

f_imp = pd.DataFrame(m['classifier'].feature_importances_, index=X.columns, columns=['score']).sort_values('score', ascending=False)
if not os.path.exists(f'{dout}/{pon_source}'):
    os.makedirs(f'{dout}/{pon_source}')

f_imp.to_csv(f'{dout}/feature_importance.csv', index=True)

## Calculate stratified performance by variant annotations 
strat_dict = {}

for anno_category in ['Func.refGene', 'ExonicFunc.refGene', 'Func.knownGene', 'ExonicFunc.knownGene', 'Func.ensGene', 'ExonicFunc.ensGene']:
    strat_pf = pd.DataFrame(
        index = meta[anno_category].unique(),
        columns = ['acc','roc','prec','recall','f1']
    )
    for anno in meta[anno_category].unique():
        var_ids = list(meta.iloc[meta[anno_category].values==anno,:]['var_id'].values)
        ysub_true = y_test.to_numpy().ravel()[[var_id in var_ids for var_id in y_test.index]]
        ysub_pred = y_test_pred[[var_id in var_ids for var_id in y_test.index]]
        ysub_prob = y_test_prob[[var_id in var_ids for var_id in y_test.index]]
        acc = metrics.balanced_accuracy_score(ysub_true, ysub_pred)
        try: # ROC cannot be calculated if only one class in y_true 
            roc = metrics.roc_auc_score(ysub_true, ysub_prob)
        except:
            roc = None
        prec = metrics.precision_score(ysub_true, ysub_pred)
        rec = metrics.recall_score(ysub_true, ysub_pred)
        f1 = metrics.f1_score(ysub_true, ysub_pred)
        strat_pf.loc[anno] = [acc, roc, prec, rec, f1]
    strat_dict[anno_category] = strat_pf

with open(f'{dout}/stratified_pf_dict.p', 'wb') as f:
    pickle.dump(strat_dict, f)

with open(f'{dout}/stratified_pf_dict.p', 'rb') as f:
    strat_dict = pickle.load(f)

annotation_type = 'Func.refGene'
somatic_rate = pd.DataFrame(index=meta[annotation_type].unique(), columns=['somatic_rate'])
for annofunc in meta[annotation_type].unique():
    meta_varids = list(meta.iloc[meta[annotation_type].values==annofunc].var_id)
    bool_in_meta = [x in meta_varids for x in y.index]
    somatic_rate.loc[annofunc,'somatic_rate'] = y.iloc[bool_in_meta,:].sum().values[0]/y.shape[0]

somatic_rate.to_csv(f'{dout}/somatic_rates_func_refgene.csv')

annotation_type = 'ExonicFunc.refGene'
somatic_rate = pd.DataFrame(index=meta[annotation_type].unique(), columns=['somatic_rate'])
for annofunc in meta[annotation_type].unique():
    meta_varids = list(meta.iloc[meta[annotation_type].values==annofunc].var_id)
    bool_in_meta = [x in meta_varids for x in y.index]
    somatic_rate.loc[annofunc,'somatic_rate'] = y.iloc[bool_in_meta,:].sum().values[0]/y.shape[0]

somatic_rate.to_csv(f'{dout}/somatic_rates_exonicfunc_refgene.csv')

## Calculate stratified performance by source: tissue-only, tissue-normal, germline 
strat_dict = {}

for source in ['tissue_normal', 'tissue_only', 'germline']:
    with open(f'{din}/{source}_var_id.txt', 'r') as f:
        var_ids = [x.strip() for x in f.readlines()]
    # 
    ysub_true = y_test.to_numpy().ravel()[[var_id in var_ids for var_id in y_test.index]]
    ysub_pred = y_test_pred[[var_id in var_ids for var_id in y_test.index]]
    ysub_prob = y_test_prob[[var_id in var_ids for var_id in y_test.index]]
    acc = metrics.balanced_accuracy_score(ysub_true, ysub_pred)
    try: # ROC cannot be calculated if only one class in y_true 
        roc = metrics.roc_auc_score(ysub_true, ysub_prob)
    except:
        roc = None
    prec = metrics.precision_score(ysub_true, ysub_pred)
    rec = metrics.recall_score(ysub_true, ysub_pred)
    f1 = metrics.f1_score(ysub_true, ysub_pred)
    strat_pf.loc[anno] = [acc, roc, prec, rec, f1]
    with open(f'{dout}/{source}_stratified_pf_dict.p', 'wb') as f:
        pickle.dump(strat_dict, f)
