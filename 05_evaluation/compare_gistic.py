import pandas as pd

regs = pd.read_csv(
    "/projects/b1131/saya/new_bbcar/data/02b_cnv/10_cleaned_cnv/reg_thres_conf90.csv",
    index_col=0,
)
amps = sum(x.startswith("Amp") for x in regs.columns)
dels = sum(x.startswith("Del") for x in regs.columns)
print(f"{amps} amplification regions and {dels} deletion regions.")

for filename in [
    "results_cyto_copy_conf90.csv",
    "results_cyto_thres_aber_conf90.csv",
    "results_reg_copy_conf90.csv",
    "results_reg_thres_conf90.csv",
]:
    data = pd.read_csv(filename, index_col=0)
    print(
        filename,
        ":",
        data.loc["Average", :]
        .iloc[[x.startswith("Test ROC") for x in data.columns]]
        .mean(),
    )
