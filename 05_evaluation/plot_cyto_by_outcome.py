import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

din = "/projects/b1131/saya/new_bbcar/data/02b_cnv"
dout = "/projects/b1131/saya/new_bbcar/plots/cnv"

cn_data = {
    "amp": pd.read_csv(
        f"{din}/10_cleaned_cnv/cyto_thres_amp.csv",
        index_col=0,
    ),
    "del": pd.read_csv(
        f"{din}/10_cleaned_cnv/cyto_thres_del.csv",
        index_col=0,
    ),
}

gene_thres = pd.read_csv(f"{din}/10_cleaned_cnv/gene_thres.csv", index_col=0)

for key, df in cn_data.items():
    cn_data[key] = (df != 0).astype(float)  # change to binary

labels = pd.read_csv("/projects/b1131/saya/new_bbcar/label_all.csv", index_col=0)
labels = labels.loc[cn_data["amp"].index, :]

amp_colors = [matplotlib.colormaps["Reds"](index) for index in np.linspace(0, 1, 3)]
del_colors = [matplotlib.colormaps["Blues"](index) for index in np.linspace(0, 1, 3)]

for cyto in ["1q21.1", "16p12.3", "6q14.1", "10q11.22"]:
    plot_data = {}

    for aber_type in ["amp", "del"]:
        data = pd.DataFrame(columns=["case", "control"], index=[0, 1])
        for intlabel, label in enumerate(["control", "case"]):
            subdata = cn_data[aber_type][cyto].loc[
                labels.iloc[labels.values == intlabel, :].index
            ]
            for aber in [0, 1]:
                data.loc[aber, label] = (
                    abs(subdata.values) == aber
                ).sum() / subdata.shape[0]
        plot_data[aber_type] = data

    fig, ax = plt.subplots(figsize=(2, 2.5))

    # Plot amplifications
    ax.bar(
        [1, 2],
        [
            plot_data["amp"]["case"].loc[1],
            plot_data["amp"]["control"].loc[1],
        ],
        edgecolor="black",
        color=amp_colors[1],
        label="Amp",
    )
    # Plot deletions
    ax.bar(
        [4, 5],
        [plot_data["del"]["case"].loc[1], plot_data["del"]["control"].loc[1]],
        edgecolor="black",
        color=del_colors[1],
        label="Del",
    )

    ax.set_xticks([1, 2, 4, 5])
    ax.set_xticklabels(["Case", "Control", "Case", "Control"], rotation=45, ha="right")
    ax.set_ylabel("Proportion")
    ax.set_ylim(0, 0.75)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    if cyto != "1q21.1":
        ax.legend()
    ax.set_title(cyto)

    plt.tight_layout()
    fig.savefig(f"{dout}/amp_del_proportions_{cyto}.png")
    plt.close()

for gene in ["BRCA1", "NBR2"]:
    plot_data = {}
    for aber_type in ["amp", "del"]:
        data = pd.DataFrame(columns=["case", "control"], index=[0, 1])
        for intlabel, label in enumerate(["control", "case"]):
            subdata = gene_thres[gene].loc[
                labels.iloc[labels.values == intlabel, :].index
            ]
            if aber_type == "amp":
                data.loc[0, label] = (subdata.values <= 0).sum() / subdata.shape[0]
                data.loc[1, label] = (subdata.values >= 1).sum() / subdata.shape[0]
            elif aber_type == "del":
                data.loc[0, label] = (subdata.values >= 0).sum() / subdata.shape[0]
                data.loc[1, label] = (subdata.values <= -1).sum() / subdata.shape[0]
        plot_data[aber_type] = data
    #
    fig, ax = plt.subplots(figsize=(2, 2.5))
    # Plot amplifications
    ax.bar(
        [1, 2],
        [
            plot_data["amp"]["case"].loc[1],
            plot_data["amp"]["control"].loc[1],
        ],
        edgecolor="black",
        color=amp_colors[1],
        label="Amp",
    )
    # Plot deletions
    ax.bar(
        [4, 5],
        [plot_data["del"]["case"].loc[1], plot_data["del"]["control"].loc[1]],
        edgecolor="black",
        color=del_colors[1],
        label="Del",
    )
    ax.set_xticks([1, 2, 4, 5])
    ax.set_xticklabels(["Case", "Control", "Case", "Control"], rotation=45, ha="right")
    ax.set_ylabel("Proportion")
    ax.set_ylim(0, 0.75)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.legend()
    ax.set_title(gene)
    plt.tight_layout()
    fig.savefig(f"{dout}/amp_del_proportions_{gene}.png")
    plt.close()
