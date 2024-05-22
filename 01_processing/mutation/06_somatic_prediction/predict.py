import os
import pickle
import sys

import numpy as np
import pandas as pd

####################
#### Load stuff ####
####################

## Data
din = sys.argv[1]
dout = sys.argv[2]

X = pd.read_csv(f"{din}/input_matched_bbcar.csv", index_col=0)
y = pd.read_csv(f"{din}/target_matched_bbcar.csv", index_col=0)

meta = pd.read_csv(
    f"{din}/04_imputed/features_imputed_bbcar.csv"
)  # this contains meta information about each variant (variant exonic function etc.)

## Model
mfn = sys.argv[3]

with open(mfn, "rb") as f:
    m = pickle.load(f)

## Indices
train_ix = pd.read_csv(f"{din}/somatic_pred_ix/bbcar/train_index.txt", header=None)
test_ix = pd.read_csv(f"{din}/somatic_pred_ix/bbcar/test_index.txt", header=None)

X_train, X_test = X.loc[train_ix.values.ravel(), :], X.loc[test_ix.values.ravel(), :]
y_train, y_test = y.loc[train_ix.values.ravel(), :], y.loc[test_ix.values.ravel(), :]

###################################
#### Predict somatic mutations ####
###################################

din = "/projects/b1131/saya/new_bbcar/data/02a_mutation/04_ml_features"
dout = "/projects/b1131/saya/new_bbcar/data/02a_mutation/07_predicted_somatic"

if not os.path.exists(f"{dout}/"):
    os.makedirs(f"{dout}/")

all_var = pd.read_csv(
    f"{din}/features_imputed.csv"
)  # this file has info on which sample has which var

X_matched = pd.read_csv(f"{din}/input_matched.csv", index_col=0)
y_matched = pd.read_csv(f"{din}/target_matched.csv", index_col=0)

X_nonmatched = pd.read_csv(f"{din}/input_nonmatched.csv", index_col=0)

print("Matrix column size matches:", str(X_matched.shape[1] == X_nonmatched.shape[1]))
print("Matrix columns matche:", str(np.all(X_matched.columns == X_nonmatched.columns)))

y_nonmatched = m.predict(X_nonmatched)

pd.DataFrame(y_nonmatched, index=X_nonmatched.index, columns=["somatic"]).to_csv(
    f"{dout}/nonmatched.csv"
)
