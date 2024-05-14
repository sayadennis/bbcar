# pylint: disable=import-error

import os

import numpy as np
import pandas as pd
from inmoose.pycombat import pycombat_norm

din = "/projects/b1131/saya/new_bbcar/data/02b_cnv/10_cleaned_cnv"
dout = "/projects/b1131/saya/new_bbcar/data/02b_cnv/11_batcheffect_removed"

if not os.path.exists(dout):
    os.makedirs(dout)

meta = pd.read_csv("/projects/b1131/saya/new_bbcar/meta.csv", index_col=0)
meta = meta.iloc[(meta.tissue.values == 1) & (meta.seq_ctrl.values != 1), :]
meta.sort_values("patient_id", inplace=True)
meta.set_index("patient_id", drop=True, inplace=True)

data = {
    "cyto_copy": pd.read_csv(f"{din}/cyto_copy_conf90.csv", index_col=0),
    "cyto_thres": pd.read_csv(f"{din}/cyto_thres_aber_conf90.csv", index_col=0),
    "reg_copy": pd.read_csv(f"{din}/reg_copy_conf90.csv", index_col=0),
    "reg_thres": pd.read_csv(f"{din}/reg_thres_conf90.csv", index_col=0),
}

meta = meta.loc[data["cyto_copy"].index, :]
batches = meta.batch.values

for key, df in data.items():
    if key.endswith("thres"):
        df = df + 1
    data_corrected = pycombat_norm(df.values.T, batches)
    data_corrected = data_corrected + (-1) * np.min(data_corrected)
    pd.DataFrame(data_corrected.T, index=df.index, columns=df.columns).to_csv(
        f"{dout}/{key}.csv"
    )
