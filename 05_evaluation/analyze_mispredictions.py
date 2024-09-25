import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

## Load predictions

din = "/projects/b1131/saya/new_bbcar/model_interpretations/breast_cancer_prediction"
dout = "/projects/b1131/saya/new_bbcar/plots"

mut_preds = pd.read_csv(f"{din}/no_fusion_mutation_orig_samplewise_pf.csv", index_col=0)
cnv_preds = pd.read_csv(f"{din}/no_fusion_cnv_cytothres_samplewise_pf.csv", index_col=0)

preds_all = pd.merge(
    mut_preds.mean(axis=1).to_frame(name="Mutation"),
    cnv_preds.mean(axis=1).to_frame(name="CNV"),
    right_index=True,
    left_index=True,
)

grouped = preds_all.groupby(["Mutation", "CNV"]).size().reset_index(name="count")
pivoted = grouped.pivot(index="Mutation", columns="CNV", values="count")

## Plot

fig, ax = plt.subplots(figsize=(2.6, 2.6))

g = sns.heatmap(
    pivoted,
    annot=True,
    fmt="d",
    linewidths=0.5,
    ax=ax,
    cbar=False,
    square=True,
)
ax.xaxis.tick_top()
ax.xaxis.set_label_position("top")

plt.tight_layout()
plt.savefig(f"{dout}/mispredictions_relplot.png")
plt.close()

## Study the samples that are difficult to get right

impossibles = list(preds_all.iloc[preds_all.sum(axis=1).values == 0, :].index)

meta = pd.read_csv("/projects/b1131/saya/new_bbcar/meta.csv")

clin = pd.read_csv(
    (
        "/projects/b1131/saya/new_bbcar/data/clinical"
        "/BBCaRDatabaseNU09B2_DATA_2023-04-04_0551.csv"
    ),
)

clin = clin[["study_id", "age", "bbxage", "benignbxyear", "followup", "hrpositive"]]
clin.index = [int(x.split("-")[0].replace("BBC", "")) for x in clin.study_id.values]

meta_tissue = meta.iloc[(meta.tissue.values == 1) & (meta.seq_ctrl.values == 0), :]
meta_tissue.set_index("patient_id", inplace=True)

merged = pd.merge(meta_tissue, clin, how="left", left_index=True, right_index=True)

print(merged.loc[impossibles, :])


# Check age and biopsy year correlation
from scipy.stats import f_oneway

# Check alignment metrics
metrics = pd.read_csv("/projects/b1131/saya/new_bbcar/metrics/alignment_metrics.csv")

for i in merged.index:
    sample_id = merged.loc[i, "sample_id"]
    merged.loc[i, "est_library_size"] = metrics.iloc[
        [x == sample_id for x in metrics.sample_id], :
    ]["ESTIMATED_LIBRARY_SIZE"].item()

for varname in ["bbxage", "benignbxyear", "est_library_size"]:
    f, p = f_oneway(
        merged.iloc[[x not in impossibles for x in merged.index], :][varname],
        merged.loc[impossibles, :][varname],
    )
    print(f"{varname}: F={f:.2f}, p={p:.3f}")
