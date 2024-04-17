"""
Processing functions to create model inputs from GISTIC2 outputs.
"""

# pylint: disable=too-many-locals
# pylint: disable=consider-using-enumerate
# pylint: disable=dangerous-default-value
# pylint: disable=too-many-nested-blocks

import glob
import os

import numpy as np
import pandas as pd

# CNV processing functions


def generate_gistic_features(gistic_dir, q_thres=0.05):
    """
    This is a function that takes the path to a GISTIC2.0 output directory
    And generates cleaned CNV features from the GISTIC outputs.
    The returned object is a set of Pandas DataFrames
    stored in a dictionary with the following keys:
        - 'gene_copy'
        - 'gene_thres'
        - 'reg_copy'
        - 'reg_thres'
        - 'cyto_copy'
        - 'cyto_thres_amp'
        - 'cyto_thres_del'
    All the Pandas DataFrames in the resulting dictionary share the following format:
        - rows : samples --indexed by sample name e.g. "1004_Tissue"
        - columns : features --indexed by feature name i.e. gene names, region names, cytoband names
    """
    ## create gene-level feature matrix
    # all_data_by_genes.txt file
    gene_copy = pd.read_csv(
        f"{gistic_dir}/all_data_by_genes.txt", sep="\t", index_col=0
    )
    gene_copy = gene_copy.T.iloc[2:, :]  # remove 'Gene ID' and 'Cytoband'
    # all_thresholded.by_genes.txt file
    gene_thres = pd.read_csv(
        f"{gistic_dir}/all_thresholded.by_genes.txt", sep="\t", index_col=0
    )
    gene_thres = gene_thres.T.iloc[2:, :]  # remove 'Gene ID' and 'Cytoband'
    ## create region-level feature matrix
    reg = pd.read_csv(f"{gistic_dir}/all_lesions.conf_90.txt", sep="\t", index_col=0)
    reg = reg[
        reg["q values"] <= q_thres
    ]  # eliminate rows with non-significant q-values
    reg_thres = reg[reg["Amplitude Threshold"] != "Actual Copy Change Given"]
    reg_copy = reg[reg["Amplitude Threshold"] == "Actual Copy Change Given"]
    reg_thres, reg_copy = reg_thres.iloc[:, 8:-1].T, reg_copy.iloc[:, 8:-1].T
    ## create cytoband-level feature matrix
    patids = reg.iloc[
        :, [all(char.isdigit() for char in colname) for colname in reg.columns]
    ].columns
    cytobands = reg["Descriptor"].unique()
    cyto_thres_amp = pd.DataFrame(None, index=patids, columns=cytobands)
    cyto_thres_del = pd.DataFrame(None, index=patids, columns=cytobands)
    cyto_copy = pd.DataFrame(None, index=patids, columns=cytobands)
    for cytoband in cytobands:
        for patid in patids:
            cn = reg.iloc[
                reg["Descriptor"].values == cytoband, [x == patid for x in reg.columns]
            ]  # CNA for that patient & cytoband
            # split the data frame slice into threshold values and copy number values
            cn_thres = cn.iloc[["CN" not in x for x in cn.index], :]
            cn_copy = cn.iloc[["CN" in x for x in cn.index], :]
            # get CN value and save it to cyto_copy
            cyto_copy.loc[patid, cytoband] = cn_copy.iloc[
                np.argmax(abs(cn_copy)), :
            ].values[0]
            # get threshold value and save it to cyto_thres
            if np.any(["Amplification" in x for x in cn_thres.index]):
                thres_amp = cn_thres.iloc[
                    ["Amplification" in x for x in cn_thres.index], :
                ]
                thres_amp = thres_amp.iloc[np.argmax(thres_amp)].values[0]
            else:
                thres_amp = 0.0
            #
            if np.any(["Deletion" in x for x in cn_thres.index]):
                thres_del = cn_thres.iloc[["Deletion" in x for x in cn_thres.index], :]
                thres_del = thres_del.iloc[np.argmax(thres_del)].values[0]
            else:
                thres_del = 0.0
            #
            cyto_thres_amp.loc[patid, cytoband] = thres_amp
            cyto_thres_del.loc[patid, cytoband] = -1 * thres_del
    for df in [cyto_copy, cyto_thres_amp, cyto_thres_del]:
        df.index = df.index.str.strip()
        df.columns = df.columns.str.strip()
    feature_dict = {
        "gene_copy": gene_copy,
        "gene_thres": gene_thres,
        "reg_copy": reg_copy,
        "reg_thres": reg_thres,
        "cyto_copy": cyto_copy,
        "cyto_thres_amp": cyto_thres_amp,
        "cyto_thres_del": cyto_thres_del,
        "cyto_thres_aber": cyto_thres_amp + abs(cyto_thres_del),
    }
    # change the indices from sample names (e.g. '1004_Tissue') to study IDs (e.g. 1004)
    for _, df in feature_dict.items():
        df.index = [int(x.split("_")[0]) for x in df.index]
        df.sort_index(axis=0, inplace=True)
    return feature_dict


def lesions_to_structured_regions(lesion_fn):
    """
    This is a function that takes a path to a GISTIC all_lesions file as its input,
    and returns a structured region table as a Pandas DataFrame object.
    The columns are as follows:
        - 'chrom' : string that looks like 'chr1' --string better than integer bc 'chrX' exists
        - 'start' : interger value of start coordinate of the region
        - 'end' : interger value of end coordinate of the region
    """
    lesions = pd.read_csv(lesion_fn, sep="\t")
    regions = pd.DataFrame(
        index=lesions.iloc[["CN" not in x for x in lesions["Unique Name"].values], :][
            "Unique Name"
        ].values,  # region names
        columns=["chrom", "start", "end"],
    )
    for region in regions.index:
        peak_limits = (
            lesions.iloc[
                lesions["Unique Name"].values == region,
                [x == "Wide Peak Limits" for x in lesions.columns],
            ]
            .values[0][0]
            .strip()
        )
        regions.loc[region, "chrom"] = peak_limits.split(":")[
            0
        ]  # get the string "chr<num>"
        regions.loc[region, "start"] = int(peak_limits.split(":")[1].split("-")[0])
        regions.loc[region, "end"] = int(
            peak_limits.split(":")[1].split("-")[1].split("(")[0]
        )
    return regions


def generate_ovelap_mx(reg1, reg2):
    """
    This is a function that takes two structured region files
    (ref. "lesions_to_structured_regions" above)
    and returns an overlap matrix with the following format:
        - rows : region names from the first region table
        - columns : region names from the second region table
        - entries : the number of bases of the overlaps, 0 if no overlap
    """
    ol_mx = pd.DataFrame(0, index=reg1.index, columns=reg2.index)
    for region in reg1.index:  # loop through reg1 regions
        chrom = reg1.loc[region, "chrom"]
        start = reg1.loc[region, "start"]
        end = reg1.loc[region, "end"]
        # generate boolean arrays that helps select relevant reg2 regions
        samechrom = reg2["chrom"].values == chrom
        start_contained = (reg2["start"] < start) & (reg2["end"] > start)
        end_contained = (reg2["start"] < end) & (reg2["end"] > end)
        # fill in the overlap matrix
        ol = []
        for i in range(len(samechrom)):  # loop through reg2 regions
            if samechrom[i]:
                if start_contained[i] & end_contained[i]:
                    ol.append(end - start)
                elif start_contained[i]:
                    ol.append(
                        (reg2["end"] - start)[i]
                    )  # looking at i'th region of the control
                elif end_contained[i]:
                    ol.append((end - reg2["start"])[i])
                else:
                    ol.append(0)
            else:
                ol.append(0)
        ol_mx.loc[region, :] = ol
    return ol_mx


def transform_gatk_to_reg_matrix(reg_table, gatk_dir, sampleset=[]):
    """
    This is a function that can transform GATK segments of given samples
    into a feature matrix that can be used for modeling.
    Descriptions of each argument:
        - reg_table (pd.DataFrame) : regions to consider w/ cols=['chrom','start','end']
        - gatk_dir (str) : path to the directory where the GATK segment files are stored
        - sampleset (list) : study IDs indicating which samples to create matrix for.
    """
    # if sampleset is None, define as all samples in the gatk_dir
    if len(sampleset) == 0:
        sampleset = np.sort([int(x.split("_")[0]) for x in os.listdir(f"{gatk_dir}")])
    # create region matrix and initialize with all zeros
    mx = pd.DataFrame(0, index=sampleset, columns=reg_table.index)
    # loop through GATK files based on gatk_dir and sampleset
    for fn in glob.glob(gatk_dir + "/*.csv"):
        studyix = int(fn.split("/")[-1].split(".")[0].split("_")[0])
        if studyix in sampleset:
            gatk = pd.read_csv(fn)
            # for each line in the GATK file
            for i in gatk.index:
                # get coordinates of the GATK segment
                chrom = gatk.loc[i, "CONTIG"]
                start = gatk.loc[i, "START"]
                end = gatk.loc[i, "END"]
                # find region with overlap
                chrom_match = np.array(
                    [x == chrom for x in reg_table["chrom"].values]
                )  # boolean as to whether chromosome matches
                start_include = (start > reg_table["start"].values) & (
                    start < reg_table["end"].values
                )  # boolean as to whether the segment start is included in the region
                end_include = (end > reg_table["start"].values) & (
                    end < reg_table["end"].values
                )  # boolean as to whether the segment start is included in the region
                if np.any((chrom_match & (start_include | end_include))):
                    match_reg = list(
                        reg_table.iloc[
                            (chrom_match & (start_include | end_include))
                        ].index
                    )
                    for reg in match_reg:
                        if abs(gatk.loc[i, "MEAN_LOG2_COPY_RATIO"]) > abs(
                            mx.loc[studyix, reg]
                        ):  # only overwrite if absolute aberrant value is bigger
                            mx.loc[studyix, reg] = (
                                gatk.loc[i, "MEAN_LOG2_COPY_RATIO"] - 1
                            )  # -1 because GISTIC and GATK uses different log2 standards
    mx.sort_index(axis=0, inplace=True)
    return mx
