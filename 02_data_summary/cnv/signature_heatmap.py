import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

sns.set_theme(color_codes=True)
from sklearn.preprocessing import StandardScaler

###################
#### Load data ####
###################

din = "/projects/b1131/saya/new_bbcar/data/02b_cnv/inhouse_signatures"
dout = "/projects/b1131/saya/new_bbcar/plots/cnv"

data = {
    "counts": pd.read_csv(
        f"{din}/seglen_ampdel_category_call_counts_per_sample.csv",
        index_col=0,
        header=0,
    ),
    "ratios": pd.read_csv(
        f"{din}/seglen_ampdel_category_call_ratios_per_sample.csv",
        index_col=0,
        header=0,
    ),
    "counts_combat": pd.read_csv(
        f"{din}/inhouse_cn_features_batcheffect_rm_combat.csv", index_col=0, header=0
    ),
}

for valtype, df in data.items():
    data[valtype] = df.drop(157, axis=0)

meta = pd.read_csv("/projects/b1131/saya/new_bbcar/meta.csv")
meta = meta.iloc[(meta.tissue.values == 1) & (meta.seq_ctrl.values == 0), :]
meta.set_index("patient_id", inplace=True)

labels = [meta.loc[pat_id, "batch"] for pat_id in data["counts"].index]

palette = sns.color_palette("colorblind")

for valtype, mx in data.items():
    scaler = StandardScaler()
    mx = pd.DataFrame(scaler.fit_transform(mx), index=mx.index, columns=mx.columns)
    sns.clustermap(
        mx,
        row_colors=meta.batch.map({1: palette[8], 2: palette[9], 3: palette[2]}),
        figsize=(13, 13),
    )
    plt.savefig(f"{dout}/clustermap_inhouse_cn_{valtype}_standardscaled.png")
    plt.close()
