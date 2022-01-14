import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

sys.path.append("bbcar_project/src")
import BBCarModelTraining

seed = 0

target = pd.read_csv("/share/fsmresfiles/bbcar/bbcar_label.csv", header=0, index_col=0)

## Xia predictive gene set features 
# Gene-level copy number values 
gene_copy = pd.read_csv("/share/fsmresfiles/bbcar/modified_features/gene_copy_conf90_xia_predictivegenes.csv", header=0, index_col=0)
X_train, X_test, y_train, y_test = train_test_split(gene_copy, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "/home/srd6051/bbcar_project/outputs/genecopy_conf90_crossval_record_xia_predictivegenes.csv")

# Gene-level threshold number values 
gene_thres = pd.read_csv("/share/fsmresfiles/bbcar/modified_features/gene_thres_conf90_xia_predictivegenes.csv", header=0, index_col=0)
X_train, X_test, y_train, y_test = train_test_split(gene_copy, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "/home/srd6051/bbcar_project/outputs/genethres_conf90_crossval_record_xia_predictivegenes.csv")


## Watkins significant regions
# Gene-level copy number values 
gene_copy = pd.read_csv("/share/fsmresfiles/bbcar/modified_features/gene_copy_conf90_watkins.csv", header=0, index_col=0)
X_train, X_test, y_train, y_test = train_test_split(gene_copy, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "/home/srd6051/bbcar_project/outputs/genecopy_conf90_crossval_record_watkins.csv")

# Gene-level threshold number values 
gene_thres = pd.read_csv("/share/fsmresfiles/bbcar/modified_features/gene_thres_conf90_watkins.csv", header=0, index_col=0)
X_train, X_test, y_train, y_test = train_test_split(gene_copy, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "/home/srd6051/bbcar_project/outputs/genethres_conf90_crossval_record_watkins.csv")


## Curtis METABRIC regions
# Gene-level copy number values 
gene_copy = pd.read_csv("/share/fsmresfiles/bbcar/modified_features/gene_copy_conf90_curtis.csv", header=0, index_col=0)
X_train, X_test, y_train, y_test = train_test_split(gene_copy, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "/home/srd6051/bbcar_project/outputs/genecopy_conf90_crossval_record_curtis.csv")

# Gene-level threshold number values 
gene_thres = pd.read_csv("/share/fsmresfiles/bbcar/modified_features/gene_thres_conf90_curtis.csv", header=0, index_col=0)
X_train, X_test, y_train, y_test = train_test_split(gene_copy, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "/home/srd6051/bbcar_project/outputs/genethres_conf90_crossval_record_curtis.csv")


## Xia used regions 
# Gene-level copy number values 
gene_copy = pd.read_csv("/share/fsmresfiles/bbcar/modified_features/gene_copy_conf90_xia.csv", header=0, index_col=0)
X_train, X_test, y_train, y_test = train_test_split(gene_copy, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "/home/srd6051/bbcar_project/outputs/genecopy_conf90_crossval_record_xia.csv")

# Gene-level threshold number values 
gene_thres = pd.read_csv("/share/fsmresfiles/bbcar/modified_features/gene_thres_conf90_xia.csv", header=0, index_col=0)
X_train, X_test, y_train, y_test = train_test_split(gene_copy, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "/home/srd6051/bbcar_project/outputs/genethres_conf90_crossval_record_xia.csv")

