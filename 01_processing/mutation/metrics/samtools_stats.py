import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

din = "/projects/b1131/saya/new_bbcar/metrics/samtools_qc"
plotdir = "/projects/b1131/saya/new_bbcar/plots"

meta = pd.read_csv("/projects/b1131/saya/new_bbcar/meta.csv")
meta = meta.iloc[meta.seq_ctrl.values == 0, :]

for status in ["untrimmed", "trimmed"]:
    meta[f"Bases mapped {status}"] = None
    meta[f"Bases duplicated {status}"] = None
    meta[f"Mean read length {status}"] = None
    meta[f"Mean read depth {status}"] = None
    for i in meta.index:
        sample_id = meta.loc[i, "sample_id"]
        metrics_fn = f"{din}_{status}/{sample_id}_metrics.txt"
        if os.path.isfile(metrics_fn):
            # Read the Samtools stats file
            with open(metrics_fn, "r", encoding="utf-8") as f:
                lines = [line.strip() for line in f.readlines()]
            # Loop through lines and collect relevant numbers
            for line in lines:
                if line.startswith("bases mapped:"):
                    bases_mapped = float(line.split(":")[-1].split()[0])
                elif line.startswith("bases duplicated:"):
                    bases_duplicated = float(line.split(":")[-1].strip())
                elif line.startswith("average length:"):
                    mean_length = float(line.split(":")[-1].strip())
            # Record metrics to dataframe
            mean_depth = (bases_mapped - bases_duplicated) / 6e7
            meta.loc[i, f"Bases mapped {status}"] = bases_mapped
            meta.loc[i, f"Bases duplicated {status}"] = bases_duplicated
            meta.loc[i, f"Mean read length {status}"] = mean_length
            meta.loc[i, f"Mean read depth {status}"] = mean_depth

low_depth = meta.iloc[meta["Mean read depth trimmed"].values < 30, :]

# Write sample IDs to potentially exclude
with open(f"{din}/exclude_samples.txt", "w") as f:
    for item in low_depth.sample_id.values:
        f.write(f"{item}\n")

#### Plot ####

# Set colors
num_points = 7
positions = np.linspace(0, 1, num_points)
viridis = plt.colormaps.get_cmap("viridis")
colors = viridis(positions)

fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(8, 8))

bins = np.arange(0, 301, 10)

for j, status in enumerate(["untrimmed", "trimmed"]):
    for batch in [1, 2, 3]:
        axs[batch - 1, j].hist(
            meta.iloc[meta.batch.values == batch, :][f"Mean read depth {status}"],
            bins=bins,
            color=colors[batch - 1],
            edgecolor="black",
        )
        axs[batch - 1, j].axvline(x=20, linestyle="dashed", color="red")
        axs[batch - 1, j].axvline(x=30, linestyle="dashed", color="orange")
        axs[batch - 1, j].set_xlabel("Mean read depth")
        axs[batch - 1, j].set_xlim(0, 310)
        if j == 0:
            axs[batch - 1, j].set_ylabel("Number of samples")

fig.text(0.03, 0.80, "Batch 1", rotation=90, va="center", ha="center", fontsize=14)
fig.text(0.03, 0.50, "Batch 2", rotation=90, va="center", ha="center", fontsize=14)
fig.text(0.03, 0.18, "Batch 3", rotation=90, va="center", ha="center", fontsize=14)

# Add column labels with plot titles
fig.text(0.34, 0.92, "Before trimming", va="center", ha="center", fontsize=14)
fig.text(0.77, 0.92, "After trimming", va="center", ha="center", fontsize=14)

fig.suptitle("Mean Read Depth Pre- and Post-trimming", fontsize=16)

plt.tight_layout()
plt.subplots_adjust(left=0.15, top=0.90)
fig.savefig(f"{plotdir}/align_mean_read_depth_histogram.png")
plt.close()
