import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

sys.path.append("bbcar_project/src")
import BBCarModelTraining

seed = 0

target = pd.read_csv("../../data/bbcar/bbcar_label.csv", header=0, index_col=0)

## Gene-level copy number values 
gene_copy = pd.read_csv("../../data/bbcar/gene_copy_conf75.csv", header=0, index_col=0)
X_train, X_test, y_train, y_test = train_test_split(gene_copy, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "~/bbcar_project/outputs/gene_copy_conf75_crossval_record_genecopy.csv")

# Gene-level threshold values
gene_thres = pd.read_csv("../../data/bbcar/gene_thres_conf75.csv", header=0, index_col=0)
X_train, X_test, y_train, y_test = train_test_split(gene_thres, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "~/bbcar_project/outputs/gene_copy_conf75_crossval_record_genethres.csv")

# Region-level copy number values 
reg_copy = pd.read_csv("../../data/bbcar/reg_copy_conf75.csv", header=0, index_col=0).T
X_train, X_test, y_train, y_test = train_test_split(reg_copy, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "~/bbcar_project/outputs/gene_copy_conf75_crossval_record_regcopy.csv")

# Region_level threshold values 
reg_thres = pd.read_csv("../../data/bbcar/reg_thres_conf75.csv", header=0, index_col=0).T
X_train, X_test, y_train, y_test = train_test_split(reg_thres, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "~/bbcar_project/outputs/gene_copy_conf75_crossval_record_regthres.csv")
