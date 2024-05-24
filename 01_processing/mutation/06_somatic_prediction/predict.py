import os
import pickle
import sys

import pandas as pd

####################
#### Load stuff ####
####################

## Data
din = sys.argv[1]
dout = sys.argv[2]

if not os.path.exists(f"{dout}/"):
    os.makedirs(f"{dout}/")

## Predictors of tissue-only variants
X_nonmatched = pd.read_csv(f"{din}/input_nonmatched.csv", index_col=0)

## Model
mfn = sys.argv[3]

with open(mfn, "rb") as f:
    m = pickle.load(f)

###################################
#### Predict somatic mutations ####
###################################

y_nonmatched = m.predict(X_nonmatched)

pd.DataFrame(y_nonmatched, index=X_nonmatched.index, columns=["somatic"]).to_csv(
    f"{dout}/nonmatched.csv"
)
