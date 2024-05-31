# pylint: disable=too-many-arguments

import itertools

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib_venn import venn2

din = "/projects/b1131/saya/new_bbcar/data/02b_cnv"
dout = "/projects/b1131/saya/new_bbcar/plots/cnv"

######################################
#### Load regions and their genes ####
######################################

regions = {
    group: pd.read_csv(
        f"{din}/09_gistic2_out_conf95_{group}/all_lesions.conf_95.txt", sep="\t"
    )
    for group in ["case", "control"]
}

for group, df in regions.items():
    regions[group] = df.iloc[
        (df["q values"].values < 0.01)
        & [x.endswith(" - CN values") for x in df["Unique Name"].values],
        :,
    ]
    regions[group]["chrom"] = [
        x.split(":")[0] for x in regions[group]["Wide Peak Limits"].values
    ]
    regions[group]["start"] = [
        int(x.split(":")[1].split("-")[0])
        for x in regions[group]["Wide Peak Limits"].values
    ]
    regions[group]["end"] = [
        int(x.split(":")[1].split("-")[1].split("(")[0])
        for x in regions[group]["Wide Peak Limits"].values
    ]

region_genes = {
    group: {
        trend: pd.read_csv(
            f"{din}/09_gistic2_out_conf95_{group}/{trend}_genes.conf_95.txt", sep="\t"
        )
        for trend in ["amp", "del"]
    }
    for group in ["case", "control"]
}

#################################################################
#### Identify overlap between case/control-recurrent regions ####
#################################################################


def check_overlap(
    chr1: str,
    start1: int,
    end1: int,
    chr2: str,
    start2: int,
    end2: int,
    min_overlap: int = 1000,
) -> bool:
    if chr1 == chr2:
        # Find the overlap boundaries
        overlap_start = max(start1, start2)
        overlap_end = min(end1, end2)
        # Calculate the overlap length
        overlap_length = overlap_end - overlap_start + 1
        # Check if the overlap length meets criteria
        overlapping = overlap_length >= min_overlap
    else:
        overlapping = False
    return overlapping


for trend in ["Amp", "Del"]:
    case_regions = regions["case"].iloc[
        [x.startswith(trend) for x in regions["case"]["Unique Name"]]
    ]
    control_regions = regions["control"].iloc[
        [x.startswith(trend) for x in regions["control"]["Unique Name"]]
    ]
    ct_overlap = 0
    for i in case_regions.index:
        for j in control_regions.index:
            if check_overlap(
                case_regions.loc[i, "chrom"],
                case_regions.loc[i, "start"],
                case_regions.loc[i, "end"],
                control_regions.loc[j, "chrom"],
                control_regions.loc[j, "start"],
                control_regions.loc[j, "end"],
            ):
                ct_overlap += 1
    case_only = case_regions.shape[0] - ct_overlap
    control_only = control_regions.shape[0] - ct_overlap
    venn = venn2(
        subsets=(case_only, control_only, ct_overlap),
        set_labels=("Case Regions", "Control Regions"),
    )
    for subset in venn.subset_labels:
        if subset:  # Check if the subset label is not None
            subset.set_fontsize(24)
    # Change font size for set labels
    for label in venn.set_labels:
        if label:  # Check if the set label is not None
            label.set_fontsize(24)
    plt.title("Amplifications" if trend == "Amp" else "Deletions", fontsize=26)
    plt.savefig(f"{dout}/case_control_region_overlap_{trend}.png")
    plt.close()

################################################################################
#### For visually identified case- or control-unique cytobands, print genes ####
################################################################################

for group, trend, cytobands in zip(
    ["case", "case", "control"],
    ["del", "amp", "amp"],
    [
        ["9q13", "10q11.23", "10q26.13", "15q13.3"],
        ["4p16.3", "2q11.2", "2p11.2", "6q14.1", "7p22.1", "14q11.2"],
        ["11p15.5"],
    ],
):
    print(f"\n\n{group} - {trend}")
    for cytoband in cytobands:
        df = region_genes[group][trend]
        aber_genes = (
            df.iloc[:, [x.startswith(cytoband) for x in df.columns]]
            .dropna()
            .iloc[4:]
            .values.ravel()
        )
        print(cytoband, aber_genes)

################################################################
#### Find case- or control-unique cytobands, identify genes ####
################################################################

genes = {
    "case": {"amp": [], "del": []},
    "control": {"amp": [], "del": []},
}

for trend in ["amp", "del"]:
    case_cytos = set(x.strip() for x in regions["case"]["Descriptor"])
    control_cytos = set(x.strip() for x in regions["control"]["Descriptor"])
    case_unique = case_cytos - control_cytos
    control_unique = control_cytos - case_cytos
    # Handle cases
    df = region_genes["case"][trend]
    aber_genes = [
        (
            df.iloc[:, [x.startswith(cytoband) for x in df.columns]]
            .dropna()
            .iloc[4:]
            .values.ravel()
        )
        for cytoband in case_unique
    ]
    aber_genes = list(itertools.chain.from_iterable(aber_genes))
    genes["case"][trend] = aber_genes
    # Handle controls
    df = region_genes["control"][trend]
    aber_genes = [
        (
            df.iloc[:, [x.startswith(cytoband) for x in df.columns]]
            .dropna()
            .iloc[4:]
            .values.ravel()
        )
        for cytoband in control_unique
    ]
    aber_genes = list(itertools.chain.from_iterable(aber_genes))
    genes["control"][trend] = aber_genes
