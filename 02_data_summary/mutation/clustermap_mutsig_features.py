import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

sns.set_theme(color_codes=True)

proj_dir = "/projects/b1131/saya/new_bbcar"
sig_results_dir = f"{proj_dir}/data/02a_mutation/08_feature_matrix"
dout = f"{proj_dir}/plots/mutation"

data = {}

for valtype in ["SBS96", "DBS78", "ID83"]:
    data[valtype] = pd.read_csv(
        f"{sig_results_dir}/raw_{valtype}_features.csv", index_col=0
    )

## Import labels
meta = pd.read_csv("/projects/b1131/saya/new_bbcar/meta.csv")
meta = meta.iloc[(meta.tissue.values == 1) & (meta.seq_ctrl.values == 0), :]
meta.set_index("patient_id", inplace=True)

impossibles = [551, 1091, 1615, 1624, 694, 1220, 710, 714, 719, 763, 1324, 1382, 1402]
meta["impossible"] = np.array([x in impossibles for x in meta.index]).astype(float)

## Plot heatmap

palette = sns.color_palette("colorblind")

for valtype, mx in data.items():
    for color_label in ["batch", "label", "impossible"]:
        scaler = StandardScaler()
        mx = pd.DataFrame(scaler.fit_transform(mx), index=mx.index, columns=mx.columns)
        # colors
        row_colors = meta.loc[mx.index, color_label].map(
            {i: palette[i - 2] for i in range(4)}
        )
        # plot
        sns.clustermap(
            mx,
            row_colors=row_colors,
            figsize=(13, 13),
            xticklabels=False,
            yticklabels=False,
        )
        plt.savefig(f"{dout}/cluster_heatmap_{valtype}_color{color_label}.png")
        plt.close()

##################
#### Plot PCA ####
##################

pca = PCA(n_components=2)
data_tf = pd.DataFrame(
    pca.fit_transform(data["SBS96"]), index=data["SBS96"].index, columns=np.arange(2)
)

fig, ax = plt.subplots()

ax.scatter(
    data_tf.iloc[[x not in impossibles for x in data_tf.index], 0],
    data_tf.iloc[[x not in impossibles for x in data_tf.index], 1],
    color="c",
)
ax.scatter(
    data_tf.iloc[[x in impossibles for x in data_tf.index], 0],
    data_tf.iloc[[x in impossibles for x in data_tf.index], 1],
    color="m",
)

plt.tight_layout()
fig.savefig("/home/srd6051/test_pca.png")
plt.close()
