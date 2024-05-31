import pickle

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import shap

din = "/projects/b1131/saya/new_bbcar/model_interpretations/breast_cancer_prediction"
dout = "/projects/b1131/saya/new_bbcar/plots"

############################################
#### Plot XGB feature importance scores ####
############################################

top_n = 5

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
    top_features = f_imp[feature].columns[
        (np.argsort(-1 * f_imp[feature].mean(axis=0)))[:top_n]
    ]
    plot_data = pd.melt(
        f_imp[feature][top_features], var_name="Feature", value_name="Importance"
    )
    # Plot importance scores
    fig, ax = plt.subplots(figsize=(4.0, 3.5))
    sns.stripplot(
        data=plot_data,
        x="Feature",
        y="Importance",
        color=colors[feature],
        ax=ax,
    )
    # ax.set_title(f"Feature Importance Scores: Top {top_n} of {feature}", fontsize=14)
    ax.set_ylim(-0.005, None)
    ax.set_xlabel("")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=30, ha="right")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.tight_layout()
    fig.savefig(f"{dout}/top_feature_importance_{feature.lower()}.png")
    plt.close()

####################
#### SHAP plots ####
####################

for category, filename, top_features in zip(
    ["cnv", "mutation"],
    ["compare_gistic/results_cyto_thres_aber_shap.p", "no_fusion_mutation_orig_shap.p"],
    [
        ["1q21.1", "16p12.3", "15q11.2", "1p36.13", "17p12", "6q14.1", "10q11.22"],
        ["T[T>C]C", "A[T>G]T", "A[T>C]C", "A[C>G]T", "G[C>G]G", "T[C>A]T", "T[C>A]C"],
    ],
):
    with open(f"{din}/{filename}", "rb") as f:
        shap_values = pickle.load(f)
    # Filter shap values by top importance features
    filtered_shap_values = shap_values[:, top_features]
    # plot
    shap.plots.beeswarm(filtered_shap_values)
    plt.tight_layout()
    plt.savefig(f"{dout}/shap_{category}_xgb.png")
    plt.close()
