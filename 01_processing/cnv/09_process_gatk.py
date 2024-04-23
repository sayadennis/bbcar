# pylint: disable=missing-module-docstring

import glob
import itertools
import sys

import pandas as pd

din = sys.argv[
    1
]  # "/projects/b1131/saya/new_bbcar/data/02b_cnv/07_called_cn_segs/tissue_only"

allmx = pd.DataFrame(
    None,
    columns=[
        "name",
        "CONTIG",
        "START",
        "END",
        "NUM_POINTS_COPY_RATIO",
        "MEAN_LOG2_COPY_RATIO",
    ],
)

for fn in glob.glob(din + "/*.called.seg"):
    sample_id = fn.split("/")[-1].split(".")[0]  # e.g. '1004'
    # Count number of rows to skip (those starting with "@")
    with open(fn, "r", encoding="utf-8") as f:
        num_skip = sum(line.startswith("@") for line in f.readlines())
    # Read segment file
    data = pd.read_csv(fn, sep="\t", skiprows=num_skip)
    # Select necessary columns, add sample name column and concatenate
    data = data.iloc[:, [colname in allmx.columns for colname in data.columns]]
    data["name"] = list(itertools.repeat(sample_id, len(data)))
    allmx = pd.concat((allmx, data))

# GISTIC and GATK uses different log2 standards - sub 1
allmx["MEAN_LOG2_COPY_RATIO"] = allmx["MEAN_LOG2_COPY_RATIO"] - 1

# Save the combined matrix
allmx.to_csv(f"{din}/gistic_input_all.tsv", sep="\t", header=False, index=False)
