import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

sys.path.append("bbcar_project/src")
import BBCarModelTraining

seed = 0

target = pd.read_csv("/share/fsmresfiles/bbcar/bbcar_label.csv", header=0, index_col=0)

## Gene-level copy number values 
gene_copy = pd.read_csv("/share/fsmresfiles/bbcar/gistic_features/gene_copy_conf90.csv", header=0, index_col=0)
X_train, X_test, y_train, y_test = train_test_split(gene_copy, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "/home/srd6051/bbcar_project/outputs/gene_copy_conf90_crossval_record.csv")

# Gene-level threshold values
gene_thres = pd.read_csv("/share/fsmresfiles/bbcar/gistic_features/gene_thres_conf90.csv", header=0, index_col=0)
X_train, X_test, y_train, y_test = train_test_split(gene_thres, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "/home/srd6051/bbcar_project/outputs/gene_thres_conf90_crossval_record.csv")

# Region-level copy number values 
reg_copy = pd.read_csv("/share/fsmresfiles/bbcar/gistic_features/reg_copy_conf90.csv", header=0, index_col=0).T
X_train, X_test, y_train, y_test = train_test_split(reg_copy, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "/home/srd6051/bbcar_project/outputs/reg_copy_conf90_crossval_record.csv")

# Region_level threshold values 
reg_thres = pd.read_csv("/share/fsmresfiles/bbcar/gistic_features/reg_thres_conf90.csv", header=0, index_col=0).T
X_train, X_test, y_train, y_test = train_test_split(reg_thres, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "/home/srd6051/bbcar_project/outputs/reg_thres_conf90_crossval_record.csv")
