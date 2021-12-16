import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

sys.path.append("bbcar_project/src")
import BBCarModelTraining

seed = 0

target = pd.read_csv("/share/fsmresfiles/bbcar/bbcar_label.csv", header=0, index_col=0)

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

## Gene-level copy number values 
gene_copy = pd.read_csv("/share/fsmresfiles/bbcar/gene_copy_conf90.csv", header=0, index_col=0)
gene_copy = select_dietlein(gene_copy)
X_train, X_test, y_train, y_test = train_test_split(gene_copy, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "/home/srd6051/bbcar_project/outputs/genecopy_conf90_crossval_record_dietlein.csv")

# Gene-level threshold values
gene_thres = pd.read_csv("/share/fsmresfiles/bbcar/gene_thres_conf90.csv", header=0, index_col=0)
X_train, X_test, y_train, y_test = train_test_split(gene_thres, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "/home/srd6051/bbcar_project/outputs/genethres_conf90_crossval_record_dietlein.csv")
