import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

sys.path.append("bbcar_project/src")
import BBCarModelTraining

seed = 0

target = pd.read_csv("../../data/bbcar/bbcar_label.csv", header=0, index_col=0)

## Gene-level copy number values 
gene_copy = pd.read_csv("../../data/bbcar/genecopy_conf90_kegg.csv", header=0, index_col=0)
X_train, X_test, y_train, y_test = train_test_split(gene_copy, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(
    X_train, y_train, X_test, y_test, 
    "~/bbcar_project/outputs/genecopy_conf90_crossval_record_kegg.csv"
)

# Gene-level threshold values
gene_thres = pd.read_csv("../../data/bbcar/genethres_conf90_kegg.csv", header=0, index_col=0)
X_train, X_test, y_train, y_test = train_test_split(gene_thres, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(
    X_train, y_train, X_test, y_test, 
    "~/bbcar_project/outputs/genethres_conf90_crossval_record_kegg.csv"
)
