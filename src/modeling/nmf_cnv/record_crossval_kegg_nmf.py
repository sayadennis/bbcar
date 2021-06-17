import os
import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

sys.path.append("bbcar_project/src")
import BBCarModelTraining
import BBCarNMF

seed = 0
k_list = [20, 40, 60, 80, 100, 150, 200]

target = pd.read_csv("../../data/bbcar/bbcar_label.csv", header=0, index_col=0)

## Gene-level threshold values 
inputfn = "../../data/bbcar/genethres_conf90_kegg.csv"
gene_thres = pd.read_csv(inputfn, header=0, index_col=0)
gene_thres = np.absolute(gene_thres) # does this do it? 
X_train, X_test, y_train, y_test = train_test_split(gene_thres, target, test_size=0.2, random_state=seed)
for k in k_list:
    F_train, F_test, _ = BBCarNMF.get_F(k, X_train, y_train, X_test, y_test)
    outfn = inputfn.split("/")[-1].split(".")[0] + "_nmf_k" + str(k) + "_crossval_record.csv"
    BBCarModelTraining.record_tuning(
        F_train, y_train, F_test, y_test, 
        os.path.join("~/bbcar_project/outputs", outfn)
    )
