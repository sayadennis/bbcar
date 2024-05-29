import pandas as pd

dout = "/projects/b1131/saya/new_bbcar/data"

## Load original feature matrices
mut_dn = "/projects/b1131/saya/new_bbcar/data/02a_mutation/08_feature_matrix"
mut_fn = "raw_SBS96_features.csv"

cnv_dn = "/projects/b1131/saya/new_bbcar/data/02b_cnv/10_cleaned_cnv"
cnv_fn = "cyto_thres_aber.csv"

mut = pd.read_csv(f"{mut_dn}/{mut_fn}", index_col=0)
cnv = pd.read_csv(f"{cnv_dn}/{cnv_fn}", index_col=0)

## Concatenate
combined = mut.join(cnv, how="inner")
combined.to_csv(f"{dout}/combined_mutation_cnv_cytothres.csv")
