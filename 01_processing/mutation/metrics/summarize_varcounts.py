# pylint: disable=line-too-long
# pylint: disable=duplicate-code

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker

var_dir = "/projects/b1131/saya/new_bbcar/data/02a_mutation/02_variant_calls"

####################
#### Get counts ####
####################

plot_data = {
    "tissue_only": {"preDPfilter": [], "postDPfilter": [], "postFFPEfilter": []},
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

##############
#### Plot ####
##############

plotlabels = {
    "preDPfilter": "Pre-filter",
    "postDPfilter": "Read depth filter",
    "postFFPEfilter": "Read depth + FFPE filters",
    "tissue_only": "Tissue Only",
    "tissue_normal": "Tissue-Normal Pairs",
    "germline_only": "Germline Only",
}


def format_with_commas(x, pos):  # pylint: disable=unused-argument
    return f"{x:,.0f}"


fig, axs = plt.subplots(nrows=1, ncols=3, figsize=(10, 6))

for i, (mode, ctdict) in enumerate(plot_data.items()):
    axs[i].violinplot(
        ctdict.values(),
        np.arange(len(ctdict)),
        showextrema=True,
        showmedians=True,
        bw_method="silverman",
    )
    axs[i].set_ylabel("Number of variants")
    axs[i].yaxis.set_major_formatter(ticker.FuncFormatter(format_with_commas))
    axs[i].set_xticks(np.arange(len(ctdict)))
    axs[i].set_xticklabels(
        [plotlabels[x] for x in ctdict.keys()], rotation=30, ha="right"
    )
    axs[i].set_title(plotlabels[mode], fontsize=14)
    axs[i].spines[["right", "top"]].set_visible(False)
    axs[i].set_ylim(0, None)

fig.suptitle("Number of variants pre- and post-filters", fontsize=16)

plt.tight_layout()
fig.savefig("/home/srd6051/varcounts_DPfilter_violin.png")
plt.close()
