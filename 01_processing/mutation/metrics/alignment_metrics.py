import pandas as pd

meta = pd.read_csv("/projects/b1131/saya/new_bbcar/meta.csv")

## Collect alignment metrics
metrics_df = pd.DataFrame(
    index=meta.index,
    columns=[
        "sample_id",
        "patient_id",
        "type",
        "batch",
        "UNPAIRED_READS_EXAMINED",
        "READ_PAIRS_EXAMINED",
        "SECONDARY_OR_SUPPLEMENTARY_RDS",
        "UNMAPPED_READS",
        "UNPAIRED_READ_DUPLICATES",
        "READ_PAIR_DUPLICATES",
        "READ_PAIR_OPTICAL_DUPLICATES",
        "PERCENT_DUPLICATION",
        "ESTIMATED_LIBRARY_SIZE",
    ],
)

for i in meta.index:
    sample_id = meta.loc[i, "sample_id"]
    patient_id = meta.loc[i, "patient_id"]
    tissue_type = "tissue" if meta.loc[i, "tissue"] == 1 else "germline"
    with open(
        (
            f"/projects/b1131/saya/new_bbcar/data/01_alignment/{tissue_type}"
            f"/metrics/{patient_id}_{tissue_type}_reads.mdup.metrics.txt"
        ),
        "r",
    ) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith("library1"):
                metrics = [float(x) for x in line.strip().split()[1:]]
                break
    metrics_df.loc[i, "sample_id"] = sample_id
    metrics_df.loc[i, "patient_id"] = patient_id
    metrics_df.loc[i, "type"] = tissue_type
    metrics_df.loc[i, "batch"] = meta.loc[i, "batch"]
    metrics_df.loc[
        i,
        [
            "UNPAIRED_READS_EXAMINED",
            "READ_PAIRS_EXAMINED",
            "SECONDARY_OR_SUPPLEMENTARY_RDS",
            "UNMAPPED_READS",
            "UNPAIRED_READ_DUPLICATES",
            "READ_PAIR_DUPLICATES",
            "READ_PAIR_OPTICAL_DUPLICATES",
            "PERCENT_DUPLICATION",
            "ESTIMATED_LIBRARY_SIZE",
        ],
    ] = metrics

metrics_df.to_csv("/home/srd6051/alignment_metrics.csv", index=False)
