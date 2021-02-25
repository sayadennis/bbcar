import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

sys.path.append("bbcar_project/src")
import BBCarModelTraining

seed = 0

target = pd.read_csv("../../data/bbcar/bbcar_label.csv", header=0, index_col=0)

## Gene-level copy number values w/ p < 0.01
gene_copy = pd.read_csv("../../data/bbcar/gene_copy_conf75_pval001.csv", header=0, index_col=0)
X_train, X_test, y_train, y_test = train_test_split(gene_copy, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "~/bbcar_project/outputs/gene_copy_conf75_crossval_record_genecopy_pval001.csv")

## Gene-level copy number values w/ p < 0.05
gene_copy = pd.read_csv("../../data/bbcar/gene_copy_conf75_pval005.csv", header=0, index_col=0)
X_train, X_test, y_train, y_test = train_test_split(gene_copy, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "~/bbcar_project/outputs/gene_copy_conf75_crossval_record_genecopy_pval005.csv")

## Gene-level copy number values w/ p < 0.10
gene_copy = pd.read_csv("../../data/bbcar/gene_copy_conf75_pval010.csv", header=0, index_col=0)
X_train, X_test, y_train, y_test = train_test_split(gene_copy, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "~/bbcar_project/outputs/gene_copy_conf75_crossval_record_genecopy_pval010.csv")

## Gene-level copy number values w/ p < 0.20
gene_copy = pd.read_csv("../../data/bbcar/gene_copy_conf75_pval020.csv", header=0, index_col=0)
X_train, X_test, y_train, y_test = train_test_split(gene_copy, target, test_size=0.2, random_state=seed)
BBCarModelTraining.record_tuning(X_train, y_train, X_test, y_test, "~/bbcar_project/outputs/gene_copy_conf75_crossval_record_genecopy_pval020.csv")

