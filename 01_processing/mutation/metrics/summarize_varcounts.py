# pylint: disable=duplicate-code
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import ticker

var_dir = "/projects/b1131/saya/new_bbcar/data/02a_mutation/02_variant_calls"
meta = pd.read_csv("/projects/b1131/saya/new_bbcar/meta.csv")

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

# Get corresponding sample IDs for each mode
mode_sample_ids = {}
for mode in ["tissue_only", "tissue_normal", "germline_only"]:
    command = f"ls {var_dir}/{mode}/*_DPfiltered.vcf"
    result = subprocess.run(
        command, shell=True, capture_output=True, text=True, check=True
    )
    files = result.stdout.splitlines()
    mode_sample_ids[mode] = [
        int(x.rsplit("/", maxsplit=1)[-1].split("_", maxsplit=1)[0]) for x in files
    ]


def format_with_commas(x, pos):  # pylint: disable=unused-argument
    return f"{x:,.0f}"


colors = [plt.get_cmap("plasma")(i / 4) for i in range(4)]

###########################################
#### Plot histograms by filtering step ####
###########################################

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

##################################################
#### Plot histograms by sample type and batch ####
##################################################

fig, axs = plt.subplots(nrows=3, ncols=4, figsize=(16, 8))

for i, (mode, ctdict) in enumerate(plot_data.items()):
    cts = (
        ctdict["predicted_somatic"]
        if mode == "tissue_only"
        else ctdict["postFFPEfilter"]
        if mode == "tissue_normal"
        else ctdict["postDPfilter"]
    )
    # plot for all
    axs[i, 0].hist(
        cts,
        bins=20,
        edgecolor="black",
        color=colors[i + 1],
    )
    axs[i, 0].set_xlabel("Number of variants")
    axs[i, 0].set_ylabel("Number of patients")
    axs[i, 0].spines[["right", "top"]].set_visible(False)
    axs[i, 0].set_xticklabels(axs[i, 0].get_xticks(), rotation=30, ha="right")
    axs[i, 0].xaxis.set_major_formatter(ticker.FuncFormatter(format_with_commas))
    # plot separately for batches
    for batch in np.arange(1, 4):
        axs[i, batch].hist(
            np.array(cts)[
                np.array(
                    [
                        x
                        in meta.iloc[meta.batch.values == (batch), :].patient_id.values
                        for x in mode_sample_ids[mode]
                    ]
                )
            ],
            bins=20,
            edgecolor="black",
            color=colors[i + 1],
        )
        axs[i, batch].set_xlabel("Number of variants")
        axs[i, batch].set_ylabel("Number of patients")
        axs[i, batch].spines[["right", "top"]].set_visible(False)
        axs[i, batch].set_xticklabels(
            axs[i, batch].get_xticks(), rotation=30, ha="right"
        )
        axs[i, batch].xaxis.set_major_formatter(
            ticker.FuncFormatter(format_with_commas)
        )

fig.suptitle("Number of variants", fontsize=16)

fig.text(0.02, 0.80, "Tissue-only", rotation=90, va="center", ha="center", fontsize=14)
fig.text(
    0.02, 0.50, "Tissue-normal", rotation=90, va="center", ha="center", fontsize=14
)
fig.text(
    0.02, 0.20, "Germline-only", rotation=90, va="center", ha="center", fontsize=14
)

fig.text(0.18, 0.925, "All", va="center", ha="center", fontsize=14)
fig.text(0.40, 0.925, "Batch 1", va="center", ha="center", fontsize=14)
fig.text(0.65, 0.925, "Batch 2", va="center", ha="center", fontsize=14)
fig.text(0.89, 0.925, "Batch 3", va="center", ha="center", fontsize=14)

plt.tight_layout()
plt.subplots_adjust(left=0.07, top=0.91)
fig.savefig("/home/srd6051/varcounts_histograms_batch.png")
plt.close()

##################################################
#### Plot histograms of just the final counts ####
##################################################

fig, axs = plt.subplots(ncols=3, figsize=(9, 3))

for i, (mode, ctdict) in enumerate(plot_data.items()):
    cts = (
        ctdict["predicted_somatic"]
        if mode == "tissue_only"
        else ctdict["postFFPEfilter"]
        if mode == "tissue_normal"
        else ctdict["postDPfilter"]
    )
    # plot for all
    axs[i].hist(
        cts,
        bins=20,
        edgecolor="black",
        color=colors[i + 1],
    )
    axs[i].set_xlabel("Number of variants")
    axs[i].set_ylabel("Number of patients")
    axs[i].spines[["right", "top"]].set_visible(False)
    axs[i].set_xticklabels(axs[i].get_xticks(), rotation=30, ha="right")
    axs[i].xaxis.set_major_formatter(ticker.FuncFormatter(format_with_commas))
    axs[i].set_title(mode.replace("_", " ").title())


plt.tight_layout()
fig.savefig("/home/srd6051/varcounts_histograms_final.png")
plt.close()
