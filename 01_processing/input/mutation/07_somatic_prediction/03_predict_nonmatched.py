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
dout = '/projects/b1131/saya/bbcar/data/02a_mutation/07_somatic_mutations'

all_var = pd.read_csv(f'{din}/04_imputed/features_imputed.csv') # this file has info on which sample has which var

X_matched = pd.read_csv(f'{din}/input_matched.csv', index_col=0)
y_matched = pd.read_csv(f'{din}/target_matched.csv', index_col=0)

X_nonmatched = pd.read_csv(f'{din}/input_nonmatched.csv', index_col=0)

## Model 
mfn = '/home/srd6051/bbcar/out/20220614_saved_best_XGB_input.p'
with open(mfn, 'rb') as f:
    m = pickle.load(f)

y_nonmatched = m.predict(X_nonmatched) # not working because avsnp column is not integers... 