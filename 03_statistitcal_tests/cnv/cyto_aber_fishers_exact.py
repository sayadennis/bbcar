import pandas as pd
from scipy.stats import fisher_exact

din = "/projects/b1131/saya/new_bbcar/data/02b_cnv"
# dout = "/projects/b1131/saya/new_bbcar/plots/cnv"

cn_data = {
    "amp": pd.read_csv(
        f"{din}/10_cleaned_cnv/cyto_thres_amp.csv",
        index_col=0,
    ),
    "del": pd.read_csv(
        f"{din}/10_cleaned_cnv/cyto_thres_del.csv",
        index_col=0,
    ),
}

gene_thres = pd.read_csv(f"{din}/10_cleaned_cnv/gene_thres.csv", index_col=0)

for key, df in cn_data.items():
    cn_data[key] = (df != 0).astype(float)  # change to binary

labels = pd.read_csv("/projects/b1131/saya/new_bbcar/label_all.csv", index_col=0)
labels = labels.loc[cn_data["amp"].index, :]

for cyto in ["1q21.1", "16p12.3", "6q14.1", "10q11.22"]:
    for aber_type in ["amp", "del"]:
        data = pd.DataFrame(columns=["case", "control"], index=[0, 1])
        for intlabel, label in enumerate(["control", "case"]):
            subdata = cn_data[aber_type][cyto].loc[
                labels.iloc[labels.values == intlabel, :].index
            ]
            for aber in [0, 1]:
                data.loc[aber, label] = (abs(subdata.values) == aber).sum()
        stat, p = fisher_exact(data)
        print(f"\n#### {cyto} -- {aber_type} ####")
        print(data)
        print(f"Fisher's Exact Test Result: stat={stat:.2f}, p={p:.1e}")


for gene in ["BRCA1", "NBR2"]:
    for aber_type in ["amp", "del"]:
        data = pd.DataFrame(columns=["case", "control"], index=[0, 1])
        for intlabel, label in enumerate(["control", "case"]):
            subdata = gene_thres[gene].loc[
                labels.iloc[labels.values == intlabel, :].index
            ]
            if aber_type == "amp":
                data.loc[0, label] = (subdata.values <= 0).sum()
                data.loc[1, label] = (subdata.values >= 1).sum()
            elif aber_type == "del":
                data.loc[0, label] = (subdata.values >= 0).sum()
                data.loc[1, label] = (subdata.values <= -1).sum()
        stat, p = fisher_exact(data)
        print(f"\n#### {gene} -- {aber_type} ####")
        print(data)
        print(f"Fisher's Exact Test Result: stat={stat:.2f}, p={p:.1e}")
