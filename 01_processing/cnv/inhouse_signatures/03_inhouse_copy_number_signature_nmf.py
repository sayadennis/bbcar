import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import MinMaxScaler

#############################
#### Load feature matrix ####
#############################

dn = "/projects/b1131/saya/bbcar/data/02b_cnv/inhouse_signatures"
plot_dn = "/projects/b1131/saya/bbcar/plots/cnv"
fn = "inhouse_cn_features_batcheffect_rm_combat.csv"

X = pd.read_csv(f"{dn}/{fn}", index_col=0)

# EXPERIMENTAL - make all non-negative
X -= np.min(X.values)

# Reorder the columns so that they make more sense
X = X[
    [
        "-:<=0.50:0-100kb",
        "-:<=0.50:100kb-1Mb",
        "-:<=0.50:1Mb-10Mb",
        "-:<=0.50:10Mb-40Mb",  # '-:<=0.50:>40Mb',
        "-:0.50-1.00:0-100kb",
        "-:0.50-1.00:100kb-1Mb",
        "-:0.50-1.00:1Mb-10Mb",
        "-:0.50-1.00:10Mb-40Mb",  # '+:0.50-1.00:>40Mb',
        "+:1.00-2.00:0-100kb",
        "+:1.00-2.00:100kb-1Mb",
        "+:1.00-2.00:1Mb-10Mb",
        "+:1.00-2.00:10Mb-40Mb",
        "+:1.00-2.00:>40Mb",
        "+:2.00-3.00:0-100kb",
        "+:2.00-3.00:100kb-1Mb",
        "+:2.00-3.00:1Mb-10Mb",
        "+:2.00-3.00:10Mb-40Mb",  # '+:2.00-3.00:>40Mb',
        "+:3.00-4.00:0-100kb",
        "+:3.00-4.00:100kb-1Mb",
        "+:3.00-4.00:1Mb-10Mb",
        "+:3.00-4.00:10Mb-40Mb",  # '+:3.00-4.00:>40Mb',
        "+:>4.00:0-100kb",
        "+:>4.00:100kb-1Mb",
        "+:>4.00:1Mb-10Mb",
        "+:>4.00:10Mb-40Mb",  # '+:>4.00:>40Mb',
    ]
]

#############################################################
#### Perform NMF and determine best number of components ####
#############################################################

k_list = np.arange(2, 16)
nmf_scores = pd.DataFrame(index=k_list, columns=["recon_error", "stability"])

for k in k_list:
    nmf = NMF(n_components=k, init="random", random_state=0, max_iter=4000)
    W = nmf.fit_transform(X)
    H = nmf.components_
    X_recon = np.dot(W, H)
    recon_error = np.linalg.norm(
        np.array(X - X_recon), ord="fro"
    )  # frobenious norm as described in the COSMIC paper
    stability = silhouette_score(
        X, np.argmax(W, axis=1)
    )  # stability is measured by silhouette score - cluster assignment determined by W matrix
    nmf_scores.loc[k, "recon_error"] = recon_error
    nmf_scores.loc[k, "stability"] = stability

scaler = MinMaxScaler()
scaler.fit(nmf_scores)
scaled_scores = pd.DataFrame(
    scaler.transform(nmf_scores), index=nmf_scores.index, columns=nmf_scores.columns
)

best_k = (
    scaled_scores["stability"] - scaled_scores["recon_error"]
).idxmax()  # index that gives the max difference

fig, ax1 = plt.subplots()

color = "tab:red"
ax1.set_xlabel("Number of components")
ax1.set_ylabel("Stability", color=color)
ax1.plot(k_list, nmf_scores["stability"], color=color, marker="s")
ax1.tick_params(axis="y", labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = "tab:blue"
ax2.set_ylabel("Reconstruction error", color=color)
ax2.plot(k_list, nmf_scores["recon_error"], color=color, marker="s")
ax2.tick_params(axis="y", labelcolor=color)

ax2.axvline(best_k, linestyle="--", color="orange")

fig.suptitle("Stability and Reconstruction Error of NMF")

fig.tight_layout()  # otherwise the right y-label is slightly clipped
fig.savefig(f"{plot_dn}/inhouse_cn_signature_stability_recon_plot.png")
plt.close()

#############################
#### Save feature matrix ####
#############################

nmf = NMF(n_components=best_k, init="random", random_state=0, max_iter=2000)
W = nmf.fit_transform(X)
H = nmf.components_

pd.DataFrame(W, index=X.index, columns=np.arange(W.shape[1])).to_csv(
    f"{dn}/inhouse_cnv_sig_per_sample.csv", header=True, index=True
)
pd.DataFrame(H, index=np.arange(H.shape[0]), columns=X.columns).T.to_csv(
    f"{dn}/inhouse_cnv_signature_features.csv", header=True, index=True
)

#####################################
#### Plot the derived signatures ####
#####################################

H = pd.DataFrame(
    H, index=[f"Signature {i+1}" for i in np.arange(H.shape[0])], columns=X.columns
)
for feature_name in [
    "-:<=0.50:>40Mb",
    "+:0.50-1.00:>40Mb",
    "+:2.00-3.00:>40Mb",
    "+:3.00-4.00:>40Mb",
    "+:>4.00:>40Mb",
]:
    if feature_name in H.columns:
        print(
            (
                f"ERROR: Feature {feature_name} is already in the H matrix. "
                "The hard-coded part of the code might be breaking."
            )
        )
    else:
        H[feature_name] = 0.0

H = H[
    [
        "-:<=0.50:0-100kb",
        "-:<=0.50:100kb-1Mb",
        "-:<=0.50:1Mb-10Mb",
        "-:<=0.50:10Mb-40Mb",
        "-:<=0.50:>40Mb",
        "-:0.50-1.00:0-100kb",
        "-:0.50-1.00:100kb-1Mb",
        "-:0.50-1.00:1Mb-10Mb",
        "-:0.50-1.00:10Mb-40Mb",
        "+:0.50-1.00:>40Mb",
        "+:1.00-2.00:0-100kb",
        "+:1.00-2.00:100kb-1Mb",
        "+:1.00-2.00:1Mb-10Mb",
        "+:1.00-2.00:10Mb-40Mb",
        "+:1.00-2.00:>40Mb",
        "+:2.00-3.00:0-100kb",
        "+:2.00-3.00:100kb-1Mb",
        "+:2.00-3.00:1Mb-10Mb",
        "+:2.00-3.00:10Mb-40Mb",
        "+:2.00-3.00:>40Mb",
        "+:3.00-4.00:0-100kb",
        "+:3.00-4.00:100kb-1Mb",
        "+:3.00-4.00:1Mb-10Mb",
        "+:3.00-4.00:10Mb-40Mb",
        "+:3.00-4.00:>40Mb",
        "+:>4.00:0-100kb",
        "+:>4.00:100kb-1Mb",
        "+:>4.00:1Mb-10Mb",
        "+:>4.00:10Mb-40Mb",
        "+:>4.00:>40Mb",
    ]
]

fig, ax = plt.subplots(H.shape[0], 1, figsize=(6, 1.5 * H.shape[0]))

colors = []
for colname in H.columns:
    if "<=0.50" in colname:
        colors.append("darkslateblue")
    elif "0.50-1.00" in colname:
        colors.append("slategray")
    elif "1.00-2.00" in colname:
        colors.append("darkseagreen")
    elif "2.00-3.00" in colname:
        colors.append("mediumorchid")
    elif "3.00-4.00" in colname:
        colors.append("gold")
    elif ">4.00" in colname:
        colors.append("crimson")

for i in range(H.shape[0]):
    ax[i].bar(np.arange(H.shape[1]), H.iloc[i, :].values.ravel(), color=colors)
    ax[i].set_yticks([0, 5, 10, 15])
    ax[i].set_ylabel(f"Sig-{i+1}")
    ax[i].set_ylim([0, 20])
    if i == (H.shape[0] - 1):
        ax[i].set_xticks(np.arange(H.shape[1]))
        ax[i].set_xticklabels(
            ["0-100kb", "100kb-1Mb", "1Mb-10Mb", "10Mb-40Mb", ">40Mb"] * 6,
            rotation=90,
            ha="center",
            fontsize=10,
        )
    else:
        ax[i].set_xticks([])
        ax[i].set_xticklabels([])


fig.suptitle("In-house CN Signatures", fontsize=14)
plt.subplots_adjust(hspace=0.05)
plt.tight_layout()
fig.savefig(f"{plot_dn}/inhouse_signature.png")
plt.close()
