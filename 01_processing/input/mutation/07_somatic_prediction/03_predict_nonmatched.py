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

## Data
din = '/projects/b1131/saya/bbcar/data/02a_mutation/04_ml_features'
dout = '/projects/b1131/saya/bbcar/data/02a_mutation/07_predicted_somatic'

all_var = pd.read_csv(f'{din}/04_imputed/features_imputed.csv') # this file has info on which sample has which var

X_matched = pd.read_csv(f'{din}/input_matched.csv', index_col=0)
y_matched = pd.read_csv(f'{din}/target_matched.csv', index_col=0)

X_nonmatched = pd.read_csv(f'{din}/input_nonmatched.csv', index_col=0)

print('Matrix column size matches:', str(X_matched.shape[1]==X_nonmatched.shape[1]))
print('Matrix columns matche:', str(np.all(X_matched.columns==X_nonmatched.columns)))

## Model 
mfn = '/projects/b1131/saya/bbcar/models/20221005_saved_best_XGB_input_matched.p'
with open(mfn, 'rb') as f:
    m = pickle.load(f)

y_nonmatched = m.predict(X_nonmatched.to_numpy())

pd.DataFrame(y_nonmatched, 
                index=X_nonmatched.index, 
                columns=['somatic']
).to_csv(f'{dout}/nonmatched.csv')
