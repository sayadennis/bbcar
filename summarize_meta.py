"""
This script is for creating the below for the overall BBCAR dataset.
    1) sample count tables to summarize the sample and patient compositions, and 
    2) argument files for Quest job arrays for downstream analysis 
"""

import os

import pandas as pd

meta_fn = "/projects/b1131/saya/new_bbcar/meta.csv"
argfiles_dir = "/projects/b1131/saya/new_bbcar/jobarray_args"
summary_dir = "/projects/b1131/saya/new_bbcar/summary"

# If directories do not exist, create
for dirname in [argfiles_dir, summary_dir]:
    if not os.path.exists(dirname):
        os.makedirs(dirname)

meta = pd.read_csv(meta_fn)
clin = pd.read_csv(
    (
        "/projects/b1131/saya/new_bbcar/data/clinical"
        "/BBCaRDatabaseNU09B2_DATA_2023-04-04_0551.csv"
    ),
)

#########################################
#### Create job array argument files ####
#########################################

patient_dict = {}

patient_dict["germline"] = list(
    meta.iloc[meta.germline.values == 1, :].patient_id.unique()
)
patient_dict["tissue"] = list(meta.iloc[meta.tissue.values == 1, :].patient_id.unique())

patient_dict["batch_1"] = list(meta.iloc[meta.batch.values == 1, :].patient_id.unique())
patient_dict["batch_2"] = list(meta.iloc[meta.batch.values == 2, :].patient_id.unique())
patient_dict["batch_3"] = list(meta.iloc[meta.batch.values == 3, :].patient_id.unique())

patient_dict["cases"] = list(meta.iloc[meta.label.values == 1, :].patient_id.unique())
patient_dict["controls"] = list(
    meta.iloc[meta.label.values == 0, :].patient_id.unique()
)

patient_dict["seq_ctrl"] = list(
    meta.iloc[meta.seq_ctrl.values == 1, :].patient_id.unique()
)

patient_dict["pon"] = list(meta.iloc[meta.pon.values == 1, :].patient_id.unique())

patient_dict["all"] = list(meta.patient_id.unique())

for groupname, patient_ids in patient_dict.items():
    print(f"{groupname}: {len(patient_ids)}")
    with open(
        f"{argfiles_dir}/patient_ids_{groupname}.txt", "w", encoding="utf-8"
    ) as f:
        for item in patient_ids:
            f.write(f"{item}\n")

#########################################
#### Create summary tables and plots ####
#########################################

summary_table = pd.DataFrame(
    0,
    index=[
        "Total patients",
        "Case patients",
        "Control patients",
        "Total samples",
        "Tissue samples",
        "Germline samples",
        "Sequencing ctrl samples",
    ],
    columns=["All", "Batch 1", "Batch 2", "Batch 3"],
)

# Fill batch specific information
for batch in [1, 2, 3]:
    summary_table.loc["Total patients", f"Batch {batch}"] = len(
        meta.iloc[meta.batch.values == batch, :].patient_id.unique()
    )
    summary_table.loc["Case patients", f"Batch {batch}"] = len(
        meta.iloc[
            (meta.batch.values == batch) & (meta.label.values == 1), :
        ].patient_id.unique()
    )
    summary_table.loc["Control patients", f"Batch {batch}"] = len(
        meta.iloc[
            (meta.batch.values == batch) & (meta.label.values == 0), :
        ].patient_id.unique()
    )
    summary_table.loc["Total samples", f"Batch {batch}"] = len(
        meta.iloc[meta.batch.values == batch, :].sample_id.unique()
    )
    summary_table.loc["Tissue samples", f"Batch {batch}"] = len(
        meta.iloc[
            (meta.batch.values == batch) & (meta.tissue.values == 1), :
        ].sample_id.unique()
    )
    summary_table.loc["Germline samples", f"Batch {batch}"] = len(
        meta.iloc[
            (meta.batch.values == batch) & (meta.germline.values == 1), :
        ].sample_id.unique()
    )
    summary_table.loc["Sequencing ctrl samples", f"Batch {batch}"] = len(
        meta.iloc[
            (meta.batch.values == batch) & (meta.seq_ctrl.values == 1), :
        ].sample_id.unique()
    )

# Fill overall information
summary_table.loc["Total patients", "All"] = len(meta.patient_id.unique())
summary_table.loc["Case patients", "All"] = len(
    meta.iloc[(meta.label.values == 1), :].patient_id.unique()
)
summary_table.loc["Control patients", "All"] = len(
    meta.iloc[(meta.label.values == 0), :].patient_id.unique()
)
summary_table.loc["Total samples", "All"] = len(meta.sample_id.unique())
summary_table.loc["Tissue samples", "All"] = len(
    meta.iloc[(meta.tissue.values == 1), :].sample_id.unique()
)
summary_table.loc["Germline samples", "All"] = len(
    meta.iloc[(meta.germline.values == 1), :].sample_id.unique()
)
summary_table.loc["Sequencing ctrl samples", "All"] = len(
    meta.iloc[(meta.seq_ctrl.values == 1), :].sample_id.unique()
)

summary_table.to_csv(
    f"{summary_dir}/patient_sample_composition.csv", index=True, header=True
)

meta = meta.iloc[meta.seq_ctrl.values == 0, :]  # Remove sequencing controls for now
slide_table = pd.DataFrame(
    index=["Case", "Control", "Total"],
    columns=["Tissue-only", "Tissue-germline", "Total"],
)
slide_table.loc["Case", "Tissue-only"] = len(
    meta.iloc[
        (meta.label.values == 1)
        & (meta.tissue.values == 1)
        & (meta.germline.values == 0),
        :,
    ].patient_id.unique()
)
slide_table.loc["Control", "Tissue-only"] = len(
    meta.iloc[
        (meta.label.values == 0)
        & (meta.tissue.values == 1)
        & (meta.germline.values == 0),
        :,
    ].patient_id.unique()
)

slide_table.loc["Case", "Tissue-germline"] = len(
    meta.iloc[
        (meta.label.values == 1)
        & (meta.tissue.values == 0)
        & (meta.germline.values == 1),
        :,
    ].patient_id.unique()
)
slide_table.loc["Control", "Tissue-germline"] = len(
    meta.iloc[
        (meta.label.values == 0)
        & (meta.tissue.values == 0)
        & (meta.germline.values == 1),
        :,
    ].patient_id.unique()
)

slide_table["Total"] = slide_table.sum(axis=1)
slide_table.loc["Total", :] = slide_table.sum(axis=0)

print(slide_table)

############################
#### Biopsy information ####
############################

clin = clin[["study_id", "age", "bbxage", "benignbxyear", "followup"]]
clin.index = [int(x.split("-")[0].replace("BBC", "")) for x in clin.study_id.values]

meta_tissue = meta.iloc[(meta.tissue.values == 1) & (meta.seq_ctrl.values == 0), :]
meta_tissue.set_index("patient_id", inplace=True)

merged = pd.merge(meta_tissue, clin, how="left", left_index=True, right_index=True)

print(
    f"Biopsy year: {int(min(merged.benignbxyear.values))}-{int(max(merged.benignbxyear))}"
)
print(f"BBX age: {merged.bbxage.mean():.1f} (SD±{merged.bbxage.std():.2f})")
