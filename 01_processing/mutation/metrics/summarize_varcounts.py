# pylint: disable=line-too-long
# pylint: disable=duplicate-code

import matplotlib.pyplot as plt
from matplotlib import ticker

var_dir = "/projects/b1131/saya/new_bbcar/data/02a_mutation/02_variant_calls"

####################
#### Get counts ####
####################

plot_data = {
    "tissue_only": {
        "preDPfilter": [],
        "postDPfilter": [],
        "postFFPEfilter": [],
        "predicted_somatic": [],
    },
    "tissue_normal": {"preDPfilter": [], "postDPfilter": [], "postFFPEfilter": []},
    "germline_only": {"preDPfilter": [], "postDPfilter": []},
}

for mode, ctdict in plot_data.items():
    for category, _ in ctdict.items():
        # Get lines of variant count metric TXT file
        with open(
            f"{var_dir}/{mode}/varcounts_{category}.txt",
            "r",
            encoding="utf-8",
        ) as f:
            lines = f.readlines()
        # Record counts to dictionary
        ctdict[category] = [int(line.strip().rsplit(maxsplit=1)[-1]) for line in lines]

#########################
#### Plot histograms ####
#########################


def format_with_commas(x, pos):  # pylint: disable=unused-argument
    return f"{x:,.0f}"


colors = [plt.get_cmap("plasma")(i / 4) for i in range(4)]

fig, axs = plt.subplots(nrows=3, ncols=4, figsize=(16, 8))

for i, (mode, ctdict) in enumerate(plot_data.items()):
    for j, (filtercat, cts) in enumerate(ctdict.items()):
        axs[i, j].hist(
            cts,
            bins=20,
            edgecolor="black",
            color=colors[i + 1],
        )
        axs[i, j].set_xlabel("Number of variants")
        axs[i, j].set_ylabel("Number of patients")
        axs[i, j].spines[["right", "top"]].set_visible(False)
        axs[i, j].set_xticklabels(axs[i, j].get_xticks(), rotation=30, ha="right")
        axs[i, j].xaxis.set_major_formatter(ticker.FuncFormatter(format_with_commas))

axs[1, 3].axis("off")
axs[2, 2].axis("off")
axs[2, 3].axis("off")

fig.suptitle("Number of variants at each step", fontsize=16)

fig.text(0.02, 0.80, "Tissue-only", rotation=90, va="center", ha="center", fontsize=14)
fig.text(
    0.02, 0.50, "Tissue-normal", rotation=90, va="center", ha="center", fontsize=14
)
fig.text(
    0.02, 0.20, "Germline-only", rotation=90, va="center", ha="center", fontsize=14
)

fig.text(0.18, 0.925, "Pre-filter", va="center", ha="center", fontsize=14)
fig.text(0.40, 0.925, "Depth filter", va="center", ha="center", fontsize=14)
fig.text(0.65, 0.925, "FFPE filter", va="center", ha="center", fontsize=14)
fig.text(0.89, 0.925, "Predicted Somatic", va="center", ha="center", fontsize=14)

plt.tight_layout()
plt.subplots_adjust(left=0.07, top=0.91)
fig.savefig("/home/srd6051/varcounts_histograms.png")
plt.close()
