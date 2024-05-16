# pylint: disable=missing-module-docstring

import os
import sys

import pandas as pd
import process_functions

din = sys.argv[1]
dout = sys.argv[2]

# din = "/projects/b1131/saya/bbcar/data/02b_cnv/09_gistic2_out_conf90"
# dout = "/projects/b1131/saya/bbcar/data/02b_cnv/10_cleaned_cnv"

cnv_dict = process_functions.generate_gistic_features(din)

if not os.path.isdir(dout):
    os.mkdir(dout)

for key, df in cnv_dict.items():
    df.to_csv(f"{dout}/{key}_conf90.csv", header=True, index=True)


###############################################################
#### Get genes that fall in amplification/deletion regions ####
###############################################################

# Load GISTIC2 outputs
amps = pd.read_csv(f"{din}/amp_genes.conf_95.txt", sep="\t", index_col=0, header=None)
dels = pd.read_csv(f"{din}/del_genes.conf_95.txt", sep="\t", index_col=0, header=None)

# Select regions where q-value < 0.05
amps = amps.iloc[
    :, [float(q) < 0.05 for q in amps.loc["q value"].values[:-1]] + [False]
]
dels = dels.iloc[
    :, [float(q) < 0.05 for q in dels.loc["q value"].values[:-1]] + [False]
]

genes = {}

for aber_type, df in zip(["Amplification", "Deletion"], [amps, dels]):
    genes[aber_type] = {}
    for col in df.columns:
        cytoband = df.loc["cytoband", col]
        peak = df.loc["wide peak boundaries", col]
        genes[aber_type][f"{cytoband} {peak}"] = list(df[col][4:].dropna().values)

top_genes = {"Amplification": [], "Deletion": []}

for aber_type in ["Amplification", "Deletion"]:
    for i, (peak_name, gene_list) in enumerate(genes[aber_type].items()):
        if i < 30:
            top_genes[aber_type] += gene_list

print(
    len(top_genes["Amplification"]),
    "amplification genes and",
    len(top_genes["Deletion"]),
    "deletion genes.",
)


with open("/home/srd6051/amp_genes.txt", "w") as f:
    for item in top_genes["Amplification"]:
        f.write(f"{item}\n")

with open("/home/srd6051/del_genes.txt", "w") as f:
    for item in top_genes["Deletion"]:
        f.write(f"{item}\n")
