# pylint: disable=wrong-import-position

import glob
import os
import sys

import pandas as pd

sys.path.append("bbcar/src/processing/filter_gistic_features")
import BBCarFilterRegions

data_dn = "/projects/b1122/saya/04_cleaned_cnv"
reg_dn = "/projects/b1122/saya/published_regions"
ref_dn = "/projects/b1122/saya/ref_gene_annotations"
out_dn = "/projects/b1122/saya/07_selected_region_cnv"

###################
#### Load data ####
###################

genecopy = pd.read_csv(os.path.join(data_dn, "gene_copy_conf90.csv"), index_col=0)
genethres = pd.read_csv(os.path.join(data_dn, "gene_thres_conf90.csv"), index_col=0)
regcopy = pd.read_csv(os.path.join(data_dn, "reg_copy_conf90.csv"), index_col=0)
regthres = pd.read_csv(os.path.join(data_dn, "reg_thres_conf90.csv"), index_col=0)

##############################################
#### Select by characterized driver genes ####
##############################################

#### KEGG ####
with open(os.path.join(reg_dn, "genes_kegg_cancer.txt"), "r", encoding="utf-8") as f:
    lines = f.readlines()

kegg_genes = []
for line in lines:
    kegg_genes.append(line.rstrip())

# filter
genecopy_kegg = genecopy.iloc[:, [x in kegg_genes for x in genecopy.columns]]
genethres_kegg = genethres.iloc[:, [x in kegg_genes for x in genethres.columns]]
print(f"Number of genes selected by KEGG: {genecopy_kegg.shape[1]}\n")
# save
genecopy_kegg.to_csv(
    os.path.join(out_dn, "gene_copy_conf90_kegg_genes.csv"), header=True, index=True
)
genethres_kegg.to_csv(
    os.path.join(out_dn, "gene_thres_conf90_kegg_genes.csv"), header=True, index=True
)

#### Dietlein cancer driver genes ####
dietlein_genes = []
with open(
    os.path.join(reg_dn, "genes_dietlein_driver.txt"), "r", encoding="utf-8"
) as f:
    lines = f.readlines()

for line in lines:
    dietlein_genes.append(line.rstrip())

# filter
genecopy_dietlein = genecopy.iloc[:, [x in dietlein_genes for x in genecopy.columns]]
genethres_dietlein = genethres.iloc[:, [x in dietlein_genes for x in genethres.columns]]
print(f"Number of genes selected by Dietlein: {genecopy_dietlein.shape[1]}\n")
# save
genecopy_dietlein.to_csv(
    os.path.join(out_dn, "gene_copy_conf90_dietlein_genes.csv"), header=True, index=True
)
genethres_dietlein.to_csv(
    os.path.join(out_dn, "gene_thres_conf90_dietlein_genes.csv"),
    header=True,
    index=True,
)

#### Xia predictive genes ####
xia_genes = []
with open(os.path.join(reg_dn, "genes_xia_predictive.txt"), "r", encoding="utf-8") as f:
    lines = f.readlines()

for line in lines:
    xia_genes.append(line.rstrip())

# filter
genecopy_xia = genecopy.iloc[:, [x in xia_genes for x in genecopy.columns]]
genethres_xia = genethres.iloc[:, [x in xia_genes for x in genethres.columns]]
print(f"Number of genes selected by Xia: {genecopy_xia.shape[1]}\n")
# save
genecopy_xia.to_csv(
    os.path.join(out_dn, "gene_copy_conf90_xia_genes.csv"), header=True, index=True
)
genethres_xia.to_csv(
    os.path.join(out_dn, "gene_thres_conf90_xia_genes.csv"), header=True, index=True
)

##############################################
#### Select by characterized risk regions ####
##############################################

genes_hg19 = os.path.join(ref_dn, "hg19_refGene.txt")
genes_hg18 = os.path.join(ref_dn, "hg18_refGene.txt")

#### Filter the gene-level features ####

for fn in glob.glob(
    os.path.join(data_dn, "gene_*_conf90.csv")
):  # gene_copy or gene_thres
    cnv_mx = pd.read_csv(fn, index_col=0)
    for regions_fn in glob.glob(os.path.join(reg_dn, "regions_*.tsv")):
        regname = regions_fn.split("_")[-1].split(".")[0]
        if regname == "bbcar":
            continue
        if regname == "watkins":
            genes_fn = genes_hg19
            filtered_cnv = BBCarFilterRegions.filter_mx(
                cnv_mx, regions_fn, genes_fn, features="genes"
            )
            filtered_cnv.to_csv(
                out_dn
                + "/"
                + fn.split("/")[-1].split(".")[0]
                + "_"
                + regions_fn.split("/")[-1].split(".")[0].split("_")[1]
                + "_ampdelselect.csv",
                header=True,
                index=True,
            )
        else:
            genes_fn = genes_hg18
            filtered_cnv = BBCarFilterRegions.filter_mx(
                cnv_mx, regions_fn, genes_fn, features="genes"
            )
            filtered_cnv.to_csv(
                out_dn
                + "/"
                + fn.split("/")[-1].split(".")[0]
                + "_"
                + regions_fn.split("/")[-1].split(".")[0].split("_")[1]
                + "_ampdelselect.csv",
                header=True,
                index=True,
            )


for fn in glob.glob(os.path.join(reg_dn, "reg_*_conf90.csv")):  # reg_copy or reg_thres
    print(f"Reading file {fn}...")
    cnv_mx = pd.read_csv(fn, index_col=0)
    cnv_mx = cnv_mx.T
    for regions_fn in glob.glob(os.path.join(reg_dn, "regions_*.tsv")):
        regname = regions_fn.split("_")[-1].split(".")[0]
        if regname == "bbcar":
            continue
        if regname == "watkins":
            genes_fn = genes_hg19
            filtered_cnv = BBCarFilterRegions.filter_mx(
                cnv_mx, regions_fn, genes_fn, features="regions"
            )
            filtered_cnv.to_csv(
                out_dn
                + "/"
                + fn.split("/")[-1].split(".")[0]
                + "_"
                + regions_fn.split("/")[-1].split(".")[0].split("_")[1]
                + "_ampdelselect.csv",
                header=True,
                index=True,
            )
        else:
            genes_fn = genes_hg18
            filtered_cnv = BBCarFilterRegions.filter_mx(
                cnv_mx, regions_fn, genes_fn, features="regions"
            )
            filtered_cnv.to_csv(
                out_dn
                + "/"
                + fn.split("/")[-1].split(".")[0]
                + "_"
                + regions_fn.split("/")[-1].split(".")[0].split("_")[1]
                + "_ampdelselect.csv",
                header=True,
                index=True,
            )
