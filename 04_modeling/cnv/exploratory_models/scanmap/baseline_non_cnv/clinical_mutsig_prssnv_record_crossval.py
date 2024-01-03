import os
import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

sys.path.append("bbcar_project/src")
import BBCarModelTraining

seed = 0

dn = "/share/fsmresfiles/bbcar"

target = pd.read_csv(os.path.join(dn, "bbcar_label_studyidindex.csv"), header=0, index_col=0)

## Clinical features
clin = pd.read_csv(os.path.join(dn, "additional_features/bbcar_clinical_info.csv"), header=0, index_col=0)
clin = clin.fillna(0)
aligned_input = clin.iloc[[x in target.index for x in clin.index],:]
aligned_target = target.iloc[[x in clin.index for x in target.index],:]
aligned_input.sort_index(inplace=True)
X_train, X_test, y_train, y_test = train_test_split(aligned_input, aligned_target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "/home/srd6051/bbcar_project/outputs/addfeat_crossval_record_clinical.csv")

# Mutational signature features (3 signatures)
mutsig = pd.read_csv(os.path.join(dn, "additional_features/bbcar_mutational_signature.tsv"), header=0, index_col=0, sep="\t")
aligned_input = mutsig.iloc[[x in target.index for x in mutsig.index],:]
aligned_target = target.iloc[[x in mutsig.index for x in target.index],:]
aligned_input.sort_index(inplace=True)
X_train, X_test, y_train, y_test = train_test_split(aligned_input, aligned_target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "/home/srd6051/bbcar_project/outputs/addfeat_crossval_record_mutsig.csv")

# 313 SNVs published based on PRS
snv = pd.read_csv(os.path.join(dn, "additional_features/bbcar_prs_snp.csv"), header=0, index_col=0)
aligned_input = snv.iloc[[x in target.index for x in snv.index],:]
aligned_input = snv.iloc[np.unique(snv.index),:]
#### find non-unique index values ####
uniquelist = np.array([])
for i in aligned_input.index:
    if i in uniquelist:
        print(i)
    else:
        uniquelist = np.append(uniquelist, i)
#### eliminate non-unique index values ####
dup_i = [np.where([x==631 for x in snv.index])[0][1], np.where([x==1389 for x in snv.index])[0][1]]
aligned_input = aligned_input.iloc[[(i not in dup_i) for i in range(aligned_input.shape[0])],:]
aligned_target = target.iloc[[x in snv.index for x in target.index],:]
aligned_input.sort_index(inplace=True)
X_train, X_test, y_train, y_test = train_test_split(aligned_input, aligned_target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "/home/srd6051/bbcar_project/outputs/addfeat_crossval_record_prssnv.csv")
