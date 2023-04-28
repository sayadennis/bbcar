import os
import sys
import glob
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
filter_type = sys.argv[2]

## Data
din = f'/projects/b1131/saya/bbcar/exploratory_winham_filter/{filter_type}/04_ml_features'
dout = f'/projects/b1131/saya/bbcar/exploratory_winham_filter/{filter_type}/07_predicted_somatic'

all_var = pd.read_csv(f'{din}/04_imputed/features_imputed_{pon_source}.csv') # this file has info on which sample has which var

X_matched = pd.read_csv(f'{din}/input_matched_{pon_source}.csv', index_col=0)
y_matched = pd.read_csv(f'{din}/target_matched_{pon_source}.csv', index_col=0)

X_nonmatched = pd.read_csv(f'{din}/input_nonmatched_{pon_source}.csv', index_col=0)

print('Matrix column size matches:', str(X_matched.shape[1]==X_nonmatched.shape[1]))
print('Matrix columns matche:', str(np.all(X_matched.columns==X_nonmatched.columns)))

## Model 
model_dir = f'/projects/b1131/saya/bbcar/model_interpretations/exploratory_winham_filter/{filter_type}'
model_filenames = [x.split('/')[-1] for x in glob.glob(f'{model_dir}/*.p')]
model_filenames.sort()
mfn = f'{model_dir}/{model_filenames[-1]}' # LATEST model (gets sorted by date)
with open(mfn, 'rb') as f:
    m = pickle.load(f)

y_nonmatched = m.predict(X_nonmatched.to_numpy())

pd.DataFrame(y_nonmatched, 
                index=X_nonmatched.index, 
                columns=['somatic']
).to_csv(f'{dout}/nonmatched_{pon_source}.csv')
