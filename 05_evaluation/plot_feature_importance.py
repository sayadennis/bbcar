import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

din = "/projects/b1131/saya/new_bbcar/model_interpretations/breast_cancer_prediction"
dout = "/projects/b1131/saya/new_bbcar/plots"

f_imp = {
    "Mutation": pd.read_csv(
        f"{din}/no_fusion_mutation_orig_feature_importances.csv", index_col=0
    ),
    "CNV": pd.read_csv(
        f"{din}/no_fusion_cnv_cytothres_feature_importances.csv", index_col=0
    ),
}

palette = sns.color_palette("colorblind")

colors = {
    "Mutation": palette[1],
    "CNV": palette[2],
}

for feature in ["Mutation", "CNV"]:
    top10_features = f_imp[feature].columns[
        (np.argsort(-1 * f_imp[feature].mean(axis=0)))[:10]
    ]
    plot_data = pd.melt(
        f_imp[feature][top10_features], var_name="Feature", value_name="Importance"
    )
    # Plot importance scores
    fig, ax = plt.subplots(figsize=(8, 4))
    sns.swarmplot(
        data=plot_data,
        x="Feature",
        y="Importance",
        color=colors[feature],
        ax=ax,
    )
    ax.set_title(f"Feature Importance Scores: Top 10 of {feature}", fontsize=14)
    ax.set_ylim(-0.005, None)
    ax.set_xlabel("")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.tight_layout()
    fig.savefig(f"{dout}/top_feature_importance_{feature.lower()}.png")
    plt.close()
