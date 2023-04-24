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
dout = '/projects/b1131/saya/bbcar/data/02a_mutation/07_predicted_somatic'

all_var = pd.read_csv(f'{din}/04_imputed/features_imputed_{pon_source}.csv') # this file has info on which sample has which var

X_matched = pd.read_csv(f'{din}/input_matched_{pon_source}.csv', index_col=0)
y_matched = pd.read_csv(f'{din}/target_matched_{pon_source}.csv', index_col=0)

X_nonmatched = pd.read_csv(f'{din}/input_nonmatched_{pon_source}.csv', index_col=0)

print('Matrix column size matches:', str(X_matched.shape[1]==X_nonmatched.shape[1]))
print('Matrix columns matche:', str(np.all(X_matched.columns==X_nonmatched.columns)))

## Model 
mfn = f'/projects/b1131/saya/bbcar/model_interpretations/{pon_source}/20230419_saved_best_XGB_input_matched_bbcar.p'
with open(mfn, 'rb') as f:
    m = pickle.load(f)

y_nonmatched = m.predict(X_nonmatched.to_numpy())

pd.DataFrame(y_nonmatched, 
                index=X_nonmatched.index, 
                columns=['somatic']
).to_csv(f'{dout}/nonmatched_{pon_source}.csv')
