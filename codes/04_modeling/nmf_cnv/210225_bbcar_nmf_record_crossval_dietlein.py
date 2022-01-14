import os
import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

sys.path.append("bbcar_project/src")
import BBCarModelTraining
import BBCarNMF

## Read the Dietlein gene set 
dietlein_geneset = []
with open("bbcar_project/dietlein_driver_genes.txt", "r") as f:
    lines = f.readlines()

for line in lines:
    dietlein_geneset.append(line.rstrip())

# Function to select only genes in the geneset
def select_dietlein(input_mx):
    dietlein_mx = input_mx.iloc[:,[x in dietlein_geneset for x in input_mx.columns]]
    return dietlein_mx


seed = 0
k_list = [20, 40, 60, 80, 100, 150, 200]

target = pd.read_csv("/share/fsmresfiles/bbcar/bbcar_label.csv", header=0, index_col=0)

## Gene-level threshold values 
inputfn = "/share/fsmresfiles/bbcar/gene_thres_conf90.csv"
gene_thres = pd.read_csv(inputfn, header=0, index_col=0)
gene_thres = np.absolute(gene_thres) # does this do it? 
gene_thres = select_dietlein(gene_thres)
X_train, X_test, y_train, y_test = train_test_split(gene_thres, target, test_size=0.2, random_state=seed)
for k in k_list:
    F_train, F_test, _ = BBCarNMF.get_F(k, X_train, y_train, X_test, y_test)
    outfn = inputfn.split("/")[-1].split(".")[0] + "_dietlein_nmf_k" + str(k) + "_crossval_record.csv"
    BBCarModelTraining.record_tuning(
        F_train, y_train, F_test, y_test, 
        os.path.join("/home/srd6051/bbcar_project/outputs", outfn)
    )
