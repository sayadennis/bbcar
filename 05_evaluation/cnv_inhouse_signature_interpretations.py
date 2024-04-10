# pylint: disable=missing-module-docstring

import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import metrics

# Define file names
dproj = "/projects/b1131/saya/bbcar"
input_dir = f"{dproj}/data/02b_cnv/inhouse_signatures"
input_fn = f"{input_dir}/inhouse_sig_batcheffect_rm_combat.csv"
target_dir = f"{dproj}/data/clinical"
target_fn = f"{target_dir}/bbcar_label_studyid_from_gatk_filenames.csv"
model_dir = (
    f"{dproj}/model_interpretations/cnv_sig_nonNMF/inhouse_sig_batcheffect_rm_combat"
)
model_fn = f"{model_dir}/20231011_saved_best_XGB_inhouse_sig_batcheffect_rm_combat.p"
indexdir = f"{dproj}/train_test_splits"

# Load data
X = pd.read_csv(input_fn, index_col=0)
y = pd.read_csv(target_fn, index_col=0)

train_ix = pd.read_csv(f"{indexdir}/train_index.txt", header=None).to_numpy().ravel()
test_ix = pd.read_csv(f"{indexdir}/test_index.txt", header=None).to_numpy().ravel()

# Split data
X_train, X_test = X.loc[train_ix, :], X.loc[test_ix, :]
y_train, y_test = y.loc[train_ix, :], y.loc[test_ix, :]

# Load trained model
with open(model_fn, "rb") as f:
    model = pickle.load(f)

# Make predictions
y_train_pred = model.predict(X_train)
y_train_prob = model.predict_proba(X_train)[:, 1]
y_test_pred = model.predict(X_test)
y_test_prob = model.predict_proba(X_test)[:, 1]


# Define simple evaluation function
def evaluate_model(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    y_prob: np.ndarray,
):
    """
    Function to obtain performance metrics from predictions.
    """
    performance = pd.DataFrame(
        index=["Overall"],
        columns=["Bal Acc", "ROC AUC", "Precision", "Recall", "F1"],
    )
    stratification = "overall"
    performance.loc[
        f"Test set {stratification}", "Bal Acc"
    ] = metrics.balanced_accuracy_score(y_true, y_pred)
    performance.loc[f"Test set {stratification}", "ROC AUC"] = metrics.roc_auc_score(
        y_true, y_prob
    )
    performance.loc[
        f"Test set {stratification}", "Precision"
    ] = metrics.precision_score(y_true, y_pred)
    performance.loc[f"Test set {stratification}", "Recall"] = metrics.recall_score(
        y_true, y_pred
    )
    performance.loc[f"Test set {stratification}", "F1"] = metrics.f1_score(
        y_true, y_pred
    )
    print(performance.T)
    return performance


results = evaluate_model(y_test.to_numpy().ravel(), y_test_pred, y_test_prob)
results = evaluate_model(y_train.to_numpy().ravel(), y_train_pred, y_train_prob)

print("Training set confusion matrix")
print(metrics.confusion_matrix(y_train, y_train_pred))

print("Test set confusion matrix")
print(metrics.confusion_matrix(y_test, y_test_pred))

# Plot feature importance scores
feature_importances = pd.DataFrame(
    model["classifier"].feature_importances_.reshape(-1, 1),
    index=X_train.columns,
    columns=["score"],
).sort_values("score", ascending=False)
feature_importances.insert(
    feature_importances.shape[1],
    "cumulative sum",
    np.cumsum(feature_importances.loc[:, "score"]),
)

top_n = 20

fig, ax = plt.subplots(2, 1, figsize=(10, 8))
ax[0].bar(
    np.arange(min(top_n, feature_importances.shape[0])),
    feature_importances.loc[:, "score"].iloc[:top_n],
)
ax[0].set_xticks(np.arange(min(top_n, feature_importances.shape[0])), [])
ax[0].set_ylabel("Feature importance")

ax[1].plot(feature_importances.loc[:, "cumulative sum"].iloc[:top_n].values)
ax[1].set_xticks(np.arange(min(top_n, feature_importances.shape[0])))
ax[1].set_xticklabels(feature_importances.index[:top_n], ha="right", rotation=30)
ax[1].set_ylabel("Cumulative feature importance")

fig.suptitle(
    f"Individual and cumulative feature importance scores (top {top_n} features)",
    size=16,
)

plt.tight_layout()
fig.savefig(
    "/projects/b1131/saya/bbcar/plots/cnv/inhouse_signature_model_feature_importance.png"
)
plt.close()
