## The purpose of this script is to modify the index of the CNV matrices and label table of BBCar data
## I am doing this so that the format is compatible with ScanMap, a supervised NMF module that Yuan developed 

import os
import sys
import numpy as np
import pandas as pd

dn = "/share/fsmresfiles/bbcar/modified_features"

orig = pd.read_csv(os.path.join(dn, "gene_thres_abs.csv"), header=0, index_col=0)
t1 = pd.read_csv(os.path.join(dn, "gene_thres_thres005_abs.csv"), header=0, index_col=0)
ts1 = pd.read_csv(os.path.join(dn, "gene_thres_thres010_abs.csv"), header=0, index_col=0)

for df in [orig, t1, ts1]:
    df["patid"] = None
    for patid in list(df.index):
        int_id = int(patid.split("_")[0])
        df.iloc[[x == patid for x in df.index],-1] = int_id
    df.set_index(["patid"], drop=True, inplace=True)
    df.sort_index(inplace=True)
    df.reset_index(inplace=True, drop=True) # reset in place


import pickle

outfn = "gene_thres_abs.pik"
with open(os.path.join(dn, outfn), "wb") as fh:
    pickle.dump(orig, fh)

outfn = "gene_thres_thres005_abs.pik"
with open(os.path.join(dn, outfn), "wb") as fh:
    pickle.dump(orig, fh)

outfn = "gene_thres_thres010_abs.pik"
with open(os.path.join(dn, outfn), "wb") as fh:
    pickle.dump(orig, fh)


## Also change index for label file 

labels = pd.read_csv("/share/fsmresfiles/bbcar/bbcar_label.csv", header=0, index_col=0)
labels["patid"] = None
for patid in list(labels.index):
    int_id = int(patid.split("_")[0])
    labels.iloc[[x == patid for x in labels.index],-1] = int_id
labels.set_index(["patid"], drop=True, inplace=True)
labels.sort_index(inplace=True)
labels.reset_index(inplace=True, drop=True)
labels.to_csv("/share/fsmresfiles/bbcar/bbcar_label_intindex.csv", header=True, index=True)

## Store train and test indices for ScanMap 
from sklearn.model_selection import train_test_split
seed = 0
df = orig # the one that has integer index 
X_train, X_test, y_train, y_test = train_test_split(df, labels, test_size=31, random_state=seed, stratify=labels)
X_train, X_val, y_train, y_val = train_test_split(X_train, y_train, test_size=30, random_state=seed, stratify=y_train)
train_index = list(X_train.index)
test_index = list(X_test.index)
val_index = list(X_val.index)

# save indices as 
with open("/share/fsmresfiles/bbcar/scanmap_data/train_indices_0.15val_0.15te.csv", "w") as fh:
    for item in train_index:
        fh.write("%s\n" % item)

with open("/share/fsmresfiles/bbcar/scanmap_data/test_indices_0.15val_0.15te.csv", "w") as fh:
    for item in test_index:
        fh.write("%s\n" % item)

with open("/share/fsmresfiles/bbcar/scanmap_data/val_indices_0.15val_0.15te.csv", "w") as fh:
    for item in val_index:
        fh.write("%s\n" % item)


## Also change index for demographic file 

clin = pd.read_csv("/share/fsmresfiles/bbcar/additional_features/bbcar_clinical_info.csv", header=0, index_col=0)
target = pd.read_csv("/share/fsmresfiles/bbcar/bbcar_label_studyidindex.csv")
for i in target["patid"].values:
    if i not in list(clin.index):
        print(i)
# this prints 467
clin_cat = pd.concat([clin, pd.DataFrame(None, index=[467], columns=clin.columns)], axis=0)
clin_new = clin_cat.iloc[[x in target["patid"].values for x in clin_cat.index],:]
clin_new.loc[467, "benign_age"] = round(np.nansum(clin_new["benign_age"].values)/200)
clin_new.loc[467, "postmeno"] = 0
clin_new.loc[467, "race_white"] = 1
clin_new.loc[467, "been_preg"] = 1
clin_new.loc[467, "rel_diag_bocancer"] = 1
clin_new.loc[467, "first_degree_rel_diag_bcancer"] = 0

clin_new.sort_index(inplace=True)
clin_new.reset_index(inplace=True, drop=True)
clin_new.to_csv("/share/fsmresfiles/bbcar/scanmap_data/pt_demo_intindex.csv", header=True, index=True)

clin = pd.read_csv("/share/fsmresfiles/bbcar/scanmap_data/pt_demo_intindex.csv", header=0, index_col=0)
mutsig = pd.read_csv("/share/fsmresfiles/bbcar/additional_features/bbcar_mutsig_4.tsv", sep="\t", header=0, index_col=0)
mutsig = mutsig.fillna(0) * 0.01
mutsig_cat = pd.concat([mutsig, pd.DataFrame(None, index=[467], columns=mutsig.columns)], axis=0)
mutsig_new = mutsig_cat.iloc[[x in labels.index for x in mutsig_cat.index],:]
mutsig_new.loc[467, "Br_J%_1"] = np.nanmean(mutsig_new["Br_J%_1"])
mutsig_new.loc[467, "Br_K%_3"] = np.nanmean(mutsig_new["Br_K%_3"])
mutsig_new.loc[467, "Br_G%_30"] = np.nanmean(mutsig_new["Br_G%_30"])
mutsig_new.loc[467, "Br_D%_MMR2"] = np.nanmean(mutsig_new["Br_D%_MMR2"])

mutsig_new.sort_index(inplace=True)
mutsig_new.reset_index(inplace=True, drop=True)

clinmut = pd.concat([clin, mutsig_new], axis=1)
clinmut.to_csv("/share/fsmresfiles/bbcar/scanmap_data/pt_demo_clinmut_intindex.csv", header=True, index=True)
