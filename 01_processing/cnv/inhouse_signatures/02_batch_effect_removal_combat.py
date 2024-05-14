# pylint: disable=import-error

import numpy as np
import pandas as pd
from inmoose.pycombat import pycombat_norm

din = "/projects/b1131/saya/new_bbcar/data/02b_cnv/inhouse_signatures"
dout = "/projects/b1131/saya/new_bbcar/data/02b_cnv/inhouse_signatures"

data = pd.read_csv(
    f"{din}/seglen_ampdel_category_call_counts_per_sample.csv", index_col=0
)

meta = pd.read_csv("/projects/b1131/saya/new_bbcar/meta.csv", index_col=0)
meta = meta.iloc[meta.tissue.values == 1, :]
batch = meta.loc[[f"{x}_tissue" for x in data.index], :].batch.values

data_corrected = pycombat_norm(data.values.T, batch)
data_corrected = data_corrected + (-1) * np.min(data_corrected)

pd.DataFrame(data_corrected.T, index=data.index, columns=data.columns).to_csv(
    f"{dout}/inhouse_cn_features_batcheffect_rm_combat.csv"
)
