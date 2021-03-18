import os
import sys
import glob
import numpy as np
import pandas as pd

sys.path.append("bbcar_project/src")
import BBCarFilterRegions

data_dn = "/path/to/data"
reg_dn = "/path/to/regions"
out_dn = "/path/to/output"

genes_hg19 = os.path.join(reg_dn, "hg19_refGene.txt")
genes_hg18 = os.path.join(reg_dn, "hg18_refGene.txt")

for fn in glob.glob(os.path.join(data_dn, "gene_*_conf90.csv")): # gene_copy or gene_thres
    cnv_mx = pd.read_csv(fn, index_col=0)
    for regions_fn in glob.glob(os.path.join(reg_dn, "regions_*.tsv")):
        regname = regions_fn.split("_")[-1].split(".")[0]
        if regname == "bbcar":
            continue
        elif regname == "watkins":
            genes_fn = genes_hg19
            filtered_cnv = BBCarFilterRegions.filter_mx(cnv_mx, regions_fn, genes_fn, features="genes")
            filtered_cnv.to_csv(
                out_dn + "/" + fn.split("/")[-1].split(".")[0] + "_" + regions_fn.split("/")[-1].split(".")[0].split("_")[1] + "_ampdelselect.csv",
                header=True, index=True
            )
        else:
            genes_fn = genes_hg18
            filtered_cnv = BBCarFilterRegions.filter_mx(cnv_mx, regions_fn, genes_fn, features="genes")
            filtered_cnv.to_csv(
                out_dn + "/" + fn.split("/")[-1].split(".")[0] + "_" + regions_fn.split("/")[-1].split(".")[0].split("_")[1] + "_ampdelselect.csv",
                header=True, index=True
            )


for fn in glob.glob(os.path.join(reg_dn, "reg_*_conf90.csv")): # reg_copy or reg_thres
    print("Reading file {}...".format(fn))
    cnv_mx = pd.read_csv(fn, index_col=0)
    cnv_mx = cnv_mx.T
    for regions_fn in glob.glob(os.path.join(reg_dn, "regions_*.tsv")):
        regname = regions_fn.split("_")[-1].split(".")[0]
        if regname == "bbcar":
            continue
        elif regname == "watkins":
            genes_fn = genes_hg19
            filtered_cnv = BBCarFilterRegions.filter_mx(cnv_mx, regions_fn, genes_fn, features="regions")
            filtered_cnv.to_csv(
                out_dn + "/" + fn.split("/")[-1].split(".")[0] + "_" + regions_fn.split("/")[-1].split(".")[0].split("_")[1] + "_ampdelselect.csv",
                header=True, index=True
            )
        else:
            genes_fn = genes_hg18
            filtered_cnv = BBCarFilterRegions.filter_mx(cnv_mx, regions_fn, genes_fn, features="regions")
            filtered_cnv.to_csv(
                out_dn + "/" + fn.split("/")[-1].split(".")[0] + "_" + regions_fn.split("/")[-1].split(".")[0].split("_")[1] + "_ampdelselect.csv",
                header=True, index=True
            )

