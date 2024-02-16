# pylint: disable=missing-module-docstring
# pylint: disable=no-member

import os
import pickle

import numpy as np
import pandas as pd
import torch
from sklearn import metrics
from suphNMF import suphNMF

###################
#### Load data ####
###################

datadir = "/projects/b1131/saya/bbcar/data"
plotdir = "/projects/b1131/saya/bbcar/plots/suphNMF_learning_curves"

X_mut = (
    pd.read_csv(
        f"{datadir}/02a_mutation/08_feature_matrix/20230423_signature_results/"
        "sbs_96_original_per_sample.csv",
        index_col=0,
    )
    * 10
)
X_cnv = (
    pd.read_csv(
        f"{datadir}/02b_cnv/inhouse_signatures/inhouse_sig_batcheffect_rm_combat.csv",
        index_col=0,
    )
    * 10
)
y = pd.read_csv(
    f"{datadir}/clinical/bbcar_label_studyid_from_gatk_filenames.csv", index_col=0
)

#####################################
#### Split out held-out test set ####
#####################################

indexdir = "/projects/b1131/saya/bbcar/train_test_splits"

train_ix = pd.read_csv(f"{indexdir}/train_index.txt", header=None).to_numpy().ravel()
test_ix = pd.read_csv(f"{indexdir}/test_index.txt", header=None).to_numpy().ravel()

X_mut_train, X_mut_test = X_mut.loc[train_ix, :], X_mut.loc[test_ix, :]
X_cnv_train, X_cnv_test = X_cnv.loc[train_ix, :], X_cnv.loc[test_ix, :]
y_train, y_test = y.loc[train_ix, :], y.loc[test_ix, :]

###########################
#### Unsupervised hNMF ####
###########################

nmf_record = {}

for k in range(4, 12):
    nmf = suphNMF(X_mut_train, X_cnv_train, y_train, n_components=k)
    cv_results, best_params = nmf.crossval_fit(
        n_iters=[1000, 2000],
        lrs=[1e-2],
        clf_weights=[0.0],
        ortho_weights=[0.0, 1e-1, 1e0, 1e1],
    )

    eval_metrics = nmf.evaluate(
        torch.tensor(X_mut_test.values, dtype=torch.float32),
        torch.tensor(X_cnv_test.values, dtype=torch.float32),
        torch.tensor(y_test.values, dtype=torch.float32),
    )

    nmf_record[k] = {
        "best_params": best_params,
        "cv_mean_roc": nmf.cv_mean_roc,
        "cv_results": cv_results,
    }

best_k = np.arange(4, 12)[
    np.argmax([nmf_record[k]["cv_mean_roc"] for k in range(4, 12)])
]
best_params = nmf_record[best_k]["best_params"]

print("\n\n\n#### Unsupervised hNMF ####")
print(f"Best k: {best_k}")
print(f"Best parameters: {best_params}")

nmf_record["best_k"] = best_k

nmf = suphNMF(X_mut_train, X_cnv_train, y_train, n_components=best_k)
nmf.fit(
    torch.tensor(X_mut_train.values, dtype=torch.float32),
    torch.tensor(X_cnv_train.values, dtype=torch.float32),
    torch.tensor(y_train.values, dtype=torch.float32),
    **best_params,
)

eval_metrics = nmf.evaluate(
    torch.tensor(X_mut_test.values, dtype=torch.float32),
    torch.tensor(X_cnv_test.values, dtype=torch.float32),
    torch.tensor(y_test.values, dtype=torch.float32),
)

y_pred = eval_metrics["y_pred"].detach().numpy()
y_score = eval_metrics["y_score"].detach().numpy()

train_roc = nmf.evaluate(
    torch.tensor(X_mut_train.values, dtype=torch.float32),
    torch.tensor(X_cnv_train.values, dtype=torch.float32),
    torch.tensor(y_train.values, dtype=torch.float32),
)["roc_auc"]
cv_mean_roc = nmf_record[best_k]["cv_mean_roc"]

print(f"\n\nTraining ROC-AUC: {train_roc}")
print(f"Cross-validation ROC-AUC: {cv_mean_roc}")
print(f"Test ROC-AUC: {metrics.roc_auc_score(y_test, y_score)}")
print(f"Test precision: {metrics.precision_score(y_test, y_pred)}")
print(f"Test recall: {metrics.recall_score(y_test, y_pred)}")
print(f"Test F1: {metrics.f1_score(y_test, y_pred)}")

dout = "/projects/b1131/saya/bbcar/model_interpretations/breast_cancer_prediction/unsupervised_hNMF"
if not os.path.exists(dout):
    os.makedirs(dout)

# Get and save W and other objects
W_test = nmf.transform(X_mut_test, X_cnv_test).detach().numpy()
W_train = nmf.transform(X_mut_train, X_cnv_train).detach().numpy()

W = pd.concat(
    (
        pd.DataFrame(W_train, index=X_mut_train.index),
        pd.DataFrame(W_test, index=X_mut_test.index),
    ),
    axis=0,
)

W.to_csv(f"{dout}/learned_W.csv", index=True, header=True)

with open(f"{dout}/nmf_record.p", "wb") as f:
    pickle.dump(nmf_record, f)

with open(f"{dout}/best_params.p", "wb") as f:
    pickle.dump(best_params, f)

with open(f"{dout}/suphNMF.p", "wb") as f:
    pickle.dump(nmf, f)

with open(f"{dout}/eval_metrics.p", "wb") as f:
    pickle.dump(eval_metrics, f)

#########################
#### Supervised hNMF ####
#########################

nmf_record = {}

for k in range(4, 12):
    nmf = suphNMF(X_mut_train, X_cnv_train, y_train, n_components=k)
    cv_results, best_params = nmf.crossval_fit(
        n_iters=[1000, 2000],
        lrs=[1e-2],
        clf_weights=[1e-1, 1e0, 1e1, 1e2],
        ortho_weights=[0.0, 1e-1, 1e0, 1e1]
        # weight_decays=[1e-3], n_iters=[1000,2000]
    )

    eval_metrics = nmf.evaluate(
        torch.tensor(X_mut_test.values, dtype=torch.float32),
        torch.tensor(X_cnv_test.values, dtype=torch.float32),
        torch.tensor(y_test.values, dtype=torch.float32),
    )

    nmf_record[k] = {
        "best_params": best_params,
        "cv_mean_roc": nmf.cv_mean_roc,
        "cv_results": cv_results,
    }

best_k = np.arange(4, 12)[
    np.argmax([nmf_record[k]["cv_mean_roc"] for k in range(4, 12)])
]
best_params = nmf_record[best_k]["best_params"]

print("\n#### Supervised hNMF ####")
print(f"Best k: {best_k}")
print(f"Best parameters: {best_params}")

nmf_record["best_k"] = best_k

nmf = suphNMF(X_mut_train, X_cnv_train, y_train, n_components=best_k)
nmf.fit(
    torch.tensor(X_mut_train.values, dtype=torch.float32),
    torch.tensor(X_cnv_train.values, dtype=torch.float32),
    torch.tensor(y_train.values, dtype=torch.float32),
    **best_params,
)

eval_metrics = nmf.evaluate(
    torch.tensor(X_mut_test.values, dtype=torch.float32),
    torch.tensor(X_cnv_test.values, dtype=torch.float32),
    torch.tensor(y_test.values, dtype=torch.float32),
)

y_pred = eval_metrics["y_pred"].detach().numpy()
y_score = eval_metrics["y_score"].detach().numpy()

train_roc = nmf.evaluate(
    torch.tensor(X_mut_train.values, dtype=torch.float32),
    torch.tensor(X_cnv_train.values, dtype=torch.float32),
    torch.tensor(y_train.values, dtype=torch.float32),
)["roc_auc"]
cv_mean_roc = nmf_record[best_k]["cv_mean_roc"]

print(f"\n\nTraining ROC-AUC: {train_roc}")
print(f"Cross-validation ROC-AUC: {cv_mean_roc}")
print(f"Test ROC-AUC: {metrics.roc_auc_score(y_test, y_score)}")
print(f"Test precision: {metrics.precision_score(y_test, y_pred)}")
print(f"Test recall: {metrics.recall_score(y_test, y_pred)}")
print(f"Test F1: {metrics.f1_score(y_test, y_pred)}")

dout = "/projects/b1131/saya/bbcar/model_interpretations/breast_cancer_prediction/supervised_hNMF"
if not os.path.exists(dout):
    os.makedirs(dout)

# Get and save W and other objects
W_test = nmf.transform(X_mut_test, X_cnv_test).detach().numpy()
W_train = nmf.transform(X_mut_train, X_cnv_train).detach().numpy()

W = pd.concat(
    (
        pd.DataFrame(W_train, index=X_mut_train.index),
        pd.DataFrame(W_test, index=X_mut_test.index),
    ),
    axis=0,
)

W.to_csv(f"{dout}/learned_W.csv", index=True, header=True)

with open(f"{dout}/nmf_record.p", "wb") as f:
    pickle.dump(nmf_record, f)

with open(f"{dout}/best_params.p", "wb") as f:
    pickle.dump(best_params, f)

with open(f"{dout}/suphNMF.p", "wb") as f:
    pickle.dump(nmf, f)

with open(f"{dout}/eval_metrics.p", "wb") as f:
    pickle.dump(eval_metrics, f)
