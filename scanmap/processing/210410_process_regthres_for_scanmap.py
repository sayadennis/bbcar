## The purpose of this script is to modify the index of the CNV matrices and label table of BBCar data
## I am doing this so that the format is compatible with ScanMap, a supervised NMF module that Yuan developed 

import os
import sys
import numpy as np
import pandas as pd
import pickle

din = "/share/fsmresfiles/bbcar/gistic_features"
dout = "/share/fsmresfiles/bbcar/modified_features"

orig = pd.read_csv(os.path.join(din, "reg_thres_conf90.csv"), header=0, index_col=0).T

#### Set row index to be integer values representing patient IDs #### 
orig["patid"] = None
for patid in list(orig.index):
    int_id = int(patid.split("_")[0])
    orig.iloc[[x == patid for x in orig.index],-1] = int_id
orig.set_index(["patid"], drop=True, inplace=True)
orig.sort_index(inplace=True)
orig.reset_index(inplace=True, drop=True) # reset in place

#### Save this into a pickle #### 
outfn = "reg_thres.pik"
with open(os.path.join(dout, outfn), "wb") as fh:
    pickle.dump(orig, fh)

