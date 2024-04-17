# pylint: disable=redefined-outer-name

import glob
import itertools
import os

import numpy as np
import pandas as pd

din = "/projects/b1131/saya/bbcar/data/02b_cnv/01_gatk_analyzed_segments"
dout = "/projects/b1131/saya/bbcar/data/02b_cnv/inhouse_signatures"

if not os.path.exists(dout):
    os.makedirs(dout)

###################
#### Load data ####
###################

data = pd.DataFrame(
    None,
    columns=[
        "name",
        "CONTIG",
        "START",
        "END",
        "NUM_POINTS_COPY_RATIO",
        "MEAN_LOG2_COPY_RATIO",
        "CALL",
    ],
)

for fn in glob.glob(din + "/*.csv"):
    samplen = fn.split("/")[-1].split(".")[0]  # example: '1004_Tissue'
    sampledata = pd.read_csv(fn)  # data = GATK segment file
    sampledata = sampledata.iloc[:, :7]  # only select columns necessary for GISTIC
    sampledata["name"] = list(
        itertools.repeat(samplen, len(sampledata))
    )  # add column indicating sample name
    data = pd.concat((data, sampledata))  # append this sample to data

data = data.iloc[data.CALL.values != "0", :]

##########################
#### Define functions ####
##########################


def seglen_bp(line: pd.Series):
    return line["END"] - line["START"]


def seg_len_category_cosmic(seglen_bp):
    if (seglen_bp > 0) & (seglen_bp <= (100 * 1e3)):
        category = "0-100kb"
    elif seglen_bp <= 1e6:
        category = "100kb-1Mb"
    elif seglen_bp <= 10 * 1e6:
        category = "1Mb-10Mb"
    elif seglen_bp <= 40 * 1e6:
        category = "10Mb-40Mb"
    else:
        category = ">40Mb"
    return category


# Add segment category by length and by copy ratio
data["SEG LENGTH BP"] = data.apply(seglen_bp, axis=1)
data["SEG LENGTH MB"] = data.apply(seglen_bp, axis=1) / 10**6
data["SEG CAT LENGTH COSMIC"] = data["SEG LENGTH BP"].apply(seg_len_category_cosmic)
seglen_bins = np.quantile(data["SEG LENGTH MB"], [0.2, 0.4, 0.6, 0.8, 0.95])
# ampdel_bins = np.quantile(2**(data.MEAN_LOG2_COPY_RATIO), [.2,.4,.6,.8,.95])
ampdel_bins = np.array([0.5, 1.0, 2.0, 3.0, 4.0])


def seg_ampdel_category(mean_log2_copy_ratio, ampdel_bins=ampdel_bins):
    larger_than = np.where([mean_log2_copy_ratio > val for val in ampdel_bins])[0]
    if len(larger_than) == len(ampdel_bins):
        lower = ampdel_bins[len(ampdel_bins) - 1]
        category = f">{lower:.2f}"
    elif len(larger_than) == 0:
        upper = ampdel_bins[0]
        category = f"<={upper:.2f}"
    else:
        lower_ix = np.max(larger_than)
        upper_ix = lower_ix + 1
        lower = ampdel_bins[lower_ix]
        upper = ampdel_bins[upper_ix]
        category = f"{lower:.2f}-{upper:.2f}"
    return category


def length_readable(mb):
    if mb < 1:
        length_string = f"{mb*(10**3):.1f}kb"
    else:
        length_string = f"{mb:.1f}Mb"
    return length_string


def seg_len_category(seglen_mb, seglen_bins=seglen_bins):
    larger_than = np.where([seglen_mb > val for val in seglen_bins])[0]
    if len(larger_than) == len(seglen_bins):  # bigger than all
        lower = seglen_bins[len(seglen_bins) - 1]
        lower = length_readable(lower)
        category = f">{lower}"
    elif len(larger_than) == 0:  # smaller than all
        upper = seglen_bins[0]
        upper = length_readable(upper)
        category = f"<={upper}"
    else:
        lower_ix = np.max(larger_than)
        upper_ix = lower_ix + 1
        lower = seglen_bins[lower_ix]
        upper = seglen_bins[upper_ix]
        lower, upper = length_readable(lower), length_readable(upper)
        category = f"{lower}-{upper}"
    return category


data["SEG CAT LENGTH ORIG"] = data["SEG LENGTH MB"].apply(seg_len_category)
data["SEG CAT LEVELS"] = (2 ** (data.MEAN_LOG2_COPY_RATIO)).apply(seg_ampdel_category)

seg_cat_name = "SEG CAT LENGTH COSMIC"

#########################################
#### Sample x length category matrix ####
#########################################

categories = data.groupby(["name", seg_cat_name, "CALL"]).size().reset_index()
counts = pd.pivot(categories, index="name", columns=(seg_cat_name, "CALL"))

# clean up
counts.index.name = None
counts.index = [x.split("_")[0] for x in counts.index]  # from '1449_Tissue' to '1449'

counts = counts.iloc[
    :, counts.columns.get_level_values("CALL") != "0"
]  # Remove CALL = 0 (no amp/del)
tuples = list(
    zip(
        counts.columns.get_level_values(seg_cat_name),
        counts.columns.get_level_values("CALL"),
    )
)
colnames = [f"{call}:{seglen}" for seglen, call in tuples]
counts.columns = colnames
counts = counts.fillna(0)
counts = counts[sorted(counts.columns)]

# calculate ratios
ratios = counts.divide(counts.sum(axis=1), axis=0)
ratios["sumcounts"] = counts.sum(axis=1)

# Save
counts.to_csv(f"{dout}/seglen_category_call_counts_per_sample.csv")
ratios.to_csv(f"{dout}/seglen_category_call_ratios_per_sample.csv")

##########################################################
#### Sample x length category x amp/del levels matrix ####
##########################################################

categories = (
    data.groupby(["name", seg_cat_name, "SEG CAT LEVELS", "CALL"]).size().reset_index()
)
counts = pd.pivot(
    categories, index="name", columns=(seg_cat_name, "SEG CAT LEVELS", "CALL")
)

# clean up
counts.index.name = None
counts.index = [x.split("_")[0] for x in counts.index]  # from '1449_Tissue' to '1449'

counts = counts.iloc[
    :, counts.columns.get_level_values("CALL") != "0"
]  # Remove CALL = 0 (no amp/del)
tuples = list(
    zip(
        counts.columns.get_level_values(seg_cat_name),
        counts.columns.get_level_values("SEG CAT LEVELS"),
        counts.columns.get_level_values("CALL"),
    )
)
colnames = [f"{call}:{seglevel}:{seglen}" for seglen, seglevel, call in tuples]
counts.columns = colnames
counts = counts.fillna(0)
counts = counts[sorted(counts.columns)]

# calculate ratios
ratios = counts.divide(counts.sum(axis=1), axis=0)
ratios["sumcounts"] = counts.sum(axis=1)

# Save
counts.to_csv(f"{dout}/seglen_ampdel_category_call_counts_per_sample.csv")
ratios.to_csv(f"{dout}/seglen_ampdel_category_call_ratios_per_sample.csv")
