# pylint: disable=missing-module-docstring

import matplotlib.pyplot as plt
import pandas as pd

###################
#### Load data ####
###################

segs = {}

tissue_types = ["Normal", "CUB", "OQ", "AN", "TU"]

for tissue_type in tissue_types:
    segs[f"methyl - {tissue_type}"] = pd.read_csv(
        f"/projects/p30791/methylation/sesame_out/copy_number/segs_{tissue_type}.tsv",
        sep="\t",
        header=None,
    )

segs["BBCAR"] = pd.read_csv(
    "/projects/b1131/saya/bbcar/data/02b_cnv/07_called_cn_segs/tissue_only/gistic_input_all.tsv",
    sep="\t",
    header=None,
)

######################################################
#### Visualize segment size and alteration levels ####
######################################################

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(9, 6))

for i, (setname, segdf) in enumerate(segs.items()):
    mean_n = segdf.groupby(0).size().mean()
    std_n = segdf.groupby(0).size().std()
    axs[i // 3, i % 3].scatter(
        (segdf[3] - segdf[2]).values / 1e6,  # x axis is segment size in Mb
        segdf[5].values,  # y axis is the CN levels
        alpha=0.15,
        s=10,
    )
    axs[i // 3, i % 3].text(
        0.25,
        0.25,
        f"{mean_n:.1f} (±{std_n:.1f})\nsegs per sample",
        transform=axs[i // 3, i % 3].transAxes,
        fontsize=9,
        verticalalignment="top",
        multialignment="center",
    )
    axs[i // 3, i % 3].spines["top"].set_visible(False)
    axs[i // 3, i % 3].spines["right"].set_visible(False)
    axs[i // 3, i % 3].set_title(setname)
    # axs[i//3, i%3].set_ylim(-4.0, 2.0)
    if i // 3 == 1:
        axs[i // 3, i % 3].set_xlabel("Segment sizes (Mb)")
    if i % 3 == 0:
        axs[i // 3, i % 3].set_ylabel("Log CN ratio estimates")

plt.tight_layout()
fig.savefig("/projects/b1131/saya/bbcar/plots/compare_seg_to_methylation.png")
plt.close()
