import pandas as pd

din = "/projects/b1131/saya/new_bbcar/data/02b_cnv"

regions = {
    group: pd.read_csv(
        f"{din}/09_gistic2_out_conf95_{group}/all_lesions.conf_95.txt", sep="\t"
    )
    for group in ["case", "control"]
}

for group, df in regions.items():
    regions[group] = df.iloc[df["q values"].values < 0.01, :]

region_genes = {
    group: {
        trend: pd.read_csv(
            f"{din}/09_gistic2_out_conf95_{group}/{trend}_genes.conf_95.txt", sep="\t"
        )
        for trend in ["amp", "del"]
    }
    for group in ["case", "control"]
}

group = "case"
trend = "del"
cytobands = ["9q13", "10q11.23", "10q26.13", "15q13.3"]
for cytoband in cytobands:
    df = region_genes[group][trend]
    aber_genes = (
        df.iloc[:, [x.startswith(cytoband) for x in df.columns]]
        .dropna()
        .iloc[4:]
        .values.ravel()
    )
    print(cytoband, aber_genes)
