import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

din = "/projects/b1131/saya/new_bbcar/metrics/samtools_qc"
plotdir = "/projects/b1131/saya/new_bbcar/plots"

meta = pd.read_csv("/projects/b1131/saya/new_bbcar/meta.csv")
meta = meta.iloc[meta.seq_ctrl.values == 0, :]

#########################################
#### Identify low mean depth samples ####
#########################################

for status in ["untrimmed", "trimmed"]:
    meta[f"Bases mapped {status}"] = None
    meta[f"Bases duplicated {status}"] = None
    meta[f"Mean read length {status}"] = None
    meta[f"Mean read depth {status}"] = None
    for i in meta.index:
        sample_id = meta.loc[i, "sample_id"]
        metrics_fn = (
            f"{din}_{status}/{sample_id}_metrics.txt"
            if status == "untrimmed"
            else f"{din}/{sample_id}_metrics.txt"
        )
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

########################################
#### Re-write job array input files ####
########################################

jobargs_dir = "/projects/b1131/saya/new_bbcar/jobarray_args"

# Read in jobarray argument files of tissue and germline samples
patient_ids = {}
for sampletype in ["tissue", "germline"]:
    with open(f"{jobargs_dir}/patient_ids_{sampletype}.txt", "r") as f:
        patient_ids[sampletype] = [line.strip() for line in f.readlines()]

print(
    "Initial patient counts: tissue",
    len(patient_ids["tissue"]),
    "// germline",
    len(patient_ids["germline"]),
)

# Remove low mean depth samples from list
for sample_id in low_depth.sample_id.values:
    patient_id, sampletype = sample_id.split("_")
    patient_ids[sampletype].remove(patient_id)

print(
    "Patient counts after removal: tissue",
    len(patient_ids["tissue"]),
    "// germline",
    len(patient_ids["germline"]),
)

# Overwrite job array argument files with the narrowed list
for sampletype, patient_id_list in patient_ids.items():
    with open(f"{jobargs_dir}/patient_ids_{sampletype}.txt", "w") as f:
        for item in patient_id_list:
            f.write(f"{item}\n")

##############
#### Plot ####
##############

# Set colors
num_points = 7
positions = np.linspace(0, 1, num_points)
viridis = plt.colormaps.get_cmap("viridis")
colors = viridis(positions)

#### Plot by batch ####

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

fig.suptitle("Mean Read Depth Pre- and Post-trimming by batch", fontsize=16)

plt.tight_layout()
plt.subplots_adjust(left=0.15, top=0.90)
fig.savefig(f"{plotdir}/align_mean_read_depth_histogram_by_batch.png")
plt.close()

#### Plot by Germline and Tissue ####

fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8, 5.3))

bins = np.arange(0, 301, 10)

for j, status in enumerate(["untrimmed", "trimmed"]):
    for i, _ in enumerate(["germline", "tissue"]):
        axs[i, j].hist(
            meta.iloc[meta.tissue.values == i, :][f"Mean read depth {status}"],
            bins=bins,
            color=colors[i + 4],
            edgecolor="black",
        )
        axs[i, j].axvline(x=20, linestyle="dashed", color="red")
        axs[i, j].axvline(x=30, linestyle="dashed", color="orange")
        axs[i, j].set_xlabel("Mean read depth")
        axs[i, j].set_xlim(0, 310)
        if j == 0:
            axs[i, j].set_ylabel("Number of samples")

fig.text(0.03, 0.66, "Germline", rotation=90, va="center", ha="center", fontsize=14)
fig.text(0.03, 0.33, "Tissue", rotation=90, va="center", ha="center", fontsize=14)

# Add column labels with plot titles
fig.text(0.34, 0.90, "Before trimming", va="center", ha="center", fontsize=14)
fig.text(0.77, 0.90, "After trimming", va="center", ha="center", fontsize=14)

fig.suptitle("Mean Read Depth Pre- and Post-trimming by tissue/germline", fontsize=16)

plt.tight_layout()
plt.subplots_adjust(left=0.15, top=0.88)
fig.savefig(f"{plotdir}/align_mean_read_depth_histogram_by_sampletype.png")
plt.close()
