import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

#######################
#### Load features ####
#######################

datadir = "/projects/b1131/saya/new_bbcar/data/02b_cnv/10_cleaned_cnv"
batchrmdir = "/projects/b1131/saya/new_bbcar/data/02b_cnv/11_batcheffect_removed"
plotdir = "/projects/b1131/saya/new_bbcar/plots/cnv"

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

# Load GISTIC2 features pre-batch effect removal
data = {
    "cyto_copy": pd.read_csv(f"{datadir}/cyto_copy_conf90.csv", index_col=0),
    "cyto_thres": pd.read_csv(f"{datadir}/cyto_thres_aber_conf90.csv", index_col=0),
    "reg_copy": pd.read_csv(f"{datadir}/reg_copy_conf90.csv", index_col=0),
    "reg_thres": pd.read_csv(f"{datadir}/reg_thres_conf90.csv", index_col=0),
}

# Load batch-effect-removed data
for key in list(data):
    data[f"{key}_batchrm"] = pd.read_csv(f"{batchrmdir}/{key}.csv", index_col=0)

meta = pd.read_csv("/projects/b1131/saya/new_bbcar/meta.csv", index_col=0)
meta = meta.iloc[(meta.tissue.values == 1) & (meta.seq_ctrl.values != 1), :]
meta.sort_values("patient_id", inplace=True)
meta.set_index("patient_id", drop=True, inplace=True)

# Exclude low-quality samples
meta = meta.loc[data["cyto_copy"].index, :]

########################################################
#### Clustermap to see patient and feature clusters ####
########################################################

palette = sns.color_palette("colorblind")

## Plot clustermaps and color rows by batch
for key, df in data.items():
    clustergrid = sns.clustermap(
        df,
        row_colors=meta.batch.map({1: palette[8], 2: palette[9], 3: palette[2]}),
    )
    plt.savefig(f"{plotdir}/clustermap_{key}.png")
    plt.close()
