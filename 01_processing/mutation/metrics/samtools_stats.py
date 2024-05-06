import os

import pandas as pd

din = "/projects/b1131/saya/new_bbcar/metrics/samtools_qc"

meta = pd.read_csv("/projects/b1131/saya/new_bbcar/meta.csv")
meta = meta.iloc[meta.seq_ctrl.values == 0, :]

meta["Bases mapped"] = None
meta["Bases duplicated"] = None
meta["Mean read length"] = None
meta["Mean read depth"] = None

for i in meta.index:
    sample_id = meta.loc[i, "sample_id"]
    metrics_fn = f"{din}/{sample_id}_metrics.txt"
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
        meta.loc[i, "Bases mapped"] = bases_mapped
        meta.loc[i, "Bases duplicated"] = bases_duplicated
        meta.loc[i, "Mean read length"] = mean_length
        meta.loc[i, "Mean read depth"] = mean_depth

low_depth = meta.iloc[meta["Mean read depth"].values < 30, :]

# Write sample IDs to potentially exclude
with open(f"{din}/exclude_samples.txt", "w") as f:
    for item in low_depth.sample_id.values:
        f.write(f"{item}\n")
