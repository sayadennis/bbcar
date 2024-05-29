import itertools

import pandas as pd

din = "/projects/b1131/saya/new_bbcar/data/02b_cnv"

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
    regions[group] = df.iloc[df["q values"].values < 0.01, :]
    regions[group] = df.iloc[
        [x.endswith(" - CN values") for x in df["Unique Name"].values], :
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
