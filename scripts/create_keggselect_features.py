import numpy as np
import pandas as pd

genecopy = pd.read_csv("../../data/bbcar/gene_copy_conf90.csv", index_col=0)
genethres = pd.read_csv("../../data/bbcar/gene_thres_conf90.csv", index_col=0)
label = pd.read_csv("../../data/bbcar/bbcar_label.csv", index_col=0)

with open("~/bbcar_project/kegg_pathways_in_cancer.txt", "r") as f:
    lines = f.readlines()

kegg_genes = []
for line in lines[2:]:
    kegg_genes.append(line.rstrip())

# note - there are 11 genes that are in the KEGG pathway that are NOT in the GISTIC features

genecopy_kegg = genecopy.iloc[:,[x in kegg_genes for x in genecopy.columns]]
genethres_kegg = genethres.iloc[:,[x in kegg_genes for x in genethres.columns]]

genecopy_kegg.to_csv("../../data/bbcar/genecopy_conf90_kegg.csv", header=True, index=True)
genethres_kegg.to_csv("../../data/bbcar/genethres_conf90_kegg.csv", header=True, index=True)
