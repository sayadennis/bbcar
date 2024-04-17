import glob
import itertools

import numpy as np
import pandas as pd

din = "/projects/b1131/saya/bbcar/data/02b_cnv/01_gatk_analyzed_segments"
dout = "/projects/b1131/saya/bbcar/data/02b_cnv/20230126_simulate_gistic_subsample/gistic_input"

for i, seed in enumerate([11, 13, 17, 19, 23]):
    np.random.seed(seed)
    flist = np.random.choice(glob.glob(din + "/*.csv"), size=26, replace=False)
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
    for fn in flist:
        samplen = fn.split("/")[-1].split(".")[0]  # example: '1004_Tissue'
        data = pd.read_csv(fn)  # data = GATK segment file
        newdata = data.iloc[:, :6]  # only select columns necessary for GISTIC
        newdata["name"] = list(
            itertools.repeat(samplen, len(newdata))
        )  # add column indicating sample name
        allmx = pd.concat((allmx, newdata))  # append this sample to allmx
    #
    allmx["MEAN_LOG2_COPY_RATIO"] = (
        allmx["MEAN_LOG2_COPY_RATIO"] - 1
    )  # GISTIC and GATK uses different log2 standards
    allmx.to_csv(
        f"{dout}/gistic_input_subset{i}.tsv", sep="\t", header=False, index=False
    )
