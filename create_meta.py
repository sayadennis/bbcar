"""
This script is for creating a new metadata file (meta.csv) 
for all the BBCAR samples, including the new ones from April 2024.
"""

from glob import glob

import pandas as pd

dout = "/projects/b1131/saya/new_bbcar"

new_meta = pd.read_excel(
    "/projects/b1131/saya/new_bbcar/data/unformatted/BBCaR_2024_sequenced_samples_correct.xlsx"
)
old_labels = pd.read_csv(
    "/projects/b1131/saya/bbcar/data/clinical/bbcar_label_studyid_from_gatk_filenames.csv",
    index_col=0,
)

##########################################
#### Create case/control label vector ####
##########################################

# Re-format patient IDs (some are string with trailing space)
new_meta["patient ID"] = [
    int(x.strip()) if isinstance(x, str) else x for x in new_meta["patient ID"].values
]

# Select only the case/control information
new_labels = new_meta[["patient ID", "Case or control"]].drop_duplicates(
    ignore_index=True
)

# Set patient ID as index and update colname
new_labels.set_index("patient ID", drop=True, inplace=True)
new_labels.index.name = None
new_labels.columns = ["label"]

# Get patient IDs of ones labeled as sequencing control
seq_controls = set(
    new_labels.iloc[new_labels.label.values == "sequencing ctrl", :].index
)

# Assess overlap
print(f"Number of sequencing controls: {len(seq_controls)}")
overlap = set(seq_controls).intersection(set(old_labels.index))
print(f"Overlap between sequencing controls and old samples: {len(overlap)}")
overlap = set(new_labels.index).intersection(set(old_labels.index))
print(f"Total overlap between old and new samples: {len(overlap)}")

# Label all samples overlapping between batch 2 and 3 as sequencing controls
for patient_id in overlap:
    new_labels.loc[patient_id, "label"] = "sequencing ctrl"

# Encode as integer
new_labels.label = new_labels.label.map({"Control": 0, "case": 1, "sequencing ctrl": 2})

# Create a label df of unique sample IDs
new_only = set(new_labels.index) - set(old_labels.index)

labels_all = pd.concat((old_labels, new_labels.loc[list(new_only), :]), axis=0)

labels_all.to_csv(f"{dout}/label_all.csv", index=True, header=True)

##############################
#### Create metadata file ####
##############################

## Get all sample and germline IDs from old batches
germline_ids = [
    int(x.split("/")[-1])
    for x in glob("/projects/b1131/saya/bbcar/data/00_raw/germline/[0-9]*")
]
tissue_ids = [
    int(x.split("/")[-1])
    for x in glob("/projects/b1131/saya/bbcar/data/00_raw/tissue/[0-9]*")
]
with open(
    "/projects/b1131/saya/bbcar/data/sample_ids_uchicago.txt", "r", encoding="utf-8"
) as f:
    uchicago_ids = [int(x.strip()) for x in f.readlines()]

meta = pd.DataFrame(
    columns=[
        "sample_id",
        "patient_id",
        "tissue",
        "germline",
        "label",
        "batch",
        "seq_ctrl",
    ],
)

## First iterate over all old samples
for patient_id in old_labels.index:
    germline = patient_id in germline_ids
    tissue = patient_id in tissue_ids
    label = old_labels.loc[patient_id, "label"]
    batch = 1 if patient_id in uchicago_ids else 2
    if tissue:
        meta = pd.concat(
            (
                meta,
                pd.DataFrame(
                    {
                        "sample_id": f"{patient_id}_tissue",
                        "patient_id": patient_id,
                        "tissue": 1,
                        "germline": 0,
                        "label": label,
                        "batch": batch,
                        "seq_ctrl": 0,  # all sequencing controls are batch 3
                    },
                    index=[0],
                ),
            )
        )
    if germline:
        meta = pd.concat(
            (
                meta,
                pd.DataFrame(
                    {
                        "sample_id": f"{patient_id}_germline",
                        "patient_id": patient_id,
                        "tissue": 0,
                        "germline": 1,
                        "label": label,
                        "batch": 2,  # NO GERMLINE ARE BATCH 1
                        "seq_ctrl": 0,  # all sequencing controls are batch 3
                    },
                    index=[0],
                ),
            )
        )


## Next iterate over all new samples
for patient_id in new_labels.index:
    germline = (
        "Saliva"
        in new_meta.iloc[new_meta["patient ID"].values == patient_id, :][
            "Tissue or Saliva"
        ].values
    )
    tissue = (
        "Tissue FFPE"
        in new_meta.iloc[new_meta["patient ID"].values == patient_id, :][
            "Tissue or Saliva"
        ].values
    )
    label = new_labels.loc[patient_id, "label"]
    # Record whether this is a sequencing control
    if label == 2:  # if sequencing control
        seq_ctrl_status = 1
        id_append = "_seqctrl"
        label = old_labels.loc[patient_id, "label"]
    else:
        seq_ctrl_status = 0
        id_append = ""
    # Record batch info (all new samples are batch 3)
    batch = 3
    if tissue:
        meta = pd.concat(
            (
                meta,
                pd.DataFrame(
                    {
                        "sample_id": f"{patient_id}_tissue{id_append}",
                        "patient_id": patient_id,
                        "tissue": 1,
                        "germline": 0,
                        "label": label,
                        "batch": 3,
                        "seq_ctrl": seq_ctrl_status,
                    },
                    index=[0],
                ),
            )
        )
    if germline:
        meta = pd.concat(
            (
                meta,
                pd.DataFrame(
                    {
                        "sample_id": f"{patient_id}_germline{id_append}",
                        "patient_id": patient_id,
                        "tissue": 0,
                        "germline": 1,
                        "label": label,
                        "batch": 3,
                        "seq_ctrl": seq_ctrl_status,
                    },
                    index=[0],
                ),
            )
        )

meta.to_csv(f"{dout}/meta.csv", index=False, header=True)
