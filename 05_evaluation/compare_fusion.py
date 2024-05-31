# pylint: disable=missing-module-docstring

from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

din = "/projects/b1131/saya/new_bbcar/model_interpretations/breast_cancer_prediction"
pf_table_fn = f"{din}/fusion_comparison.csv"

#####################################
#### Load performance dataframes ####
#####################################

performances = {
    "No fusion mutation": pd.read_csv(
        f"{din}/no_fusion_mutation_orig.csv", index_col=0
    ),
    "No fusion cnv": pd.read_csv(
        f"{din}/no_fusion_cnv_inhouse_cn_features.csv", index_col=0
    ),
    "Early fusion": pd.read_csv(f"{din}/early_fusion_cnv_inhouse.csv", index_col=0),
    "Late fusion": pd.read_csv(f"{din}/late_fusion_cnv_inhouse.csv", index_col=0),
    "Mid fusion supervised": pd.read_csv(
        f"{din}/mid_fusion_supervised_cnv_cytothres.csv",
        index_col=0,
    ),
    # "Mid fusion unsupervised": pd.read_csv(
    #    f"{din}/mid_fusion_unsupervised.csv", index_col=0
    # ),
}

## Create and print the performance comparison table for slides
compare_fusion = pd.DataFrame(
    index=performances.keys(),
    columns=["CV ROC", "Test ROC", "Test Precision", "Test Recall", "Test F1"],
)

for key, df in performances.items():
    # Get CV ROC
    mean_score = np.mean(
        df.loc["Average", [x.startswith("CV ROC") for x in df.columns]]
    )
    std_score = np.std(df.loc["Average", [x.startswith("CV ROC") for x in df.columns]])
    compare_fusion.loc[key, "CV ROC"] = f"{mean_score:.3f} (±{std_score:.3f})"
    # Get Test metrics
    for metric in ["ROC", "Precision", "Recall", "F1"]:
        mean_score = np.mean(
            df.loc["Average", [x.startswith(f"Test {metric}") for x in df.columns]]
        )
        std_score = np.std(
            df.loc["Average", [x.startswith(f"Test {metric}") for x in df.columns]]
        )
        compare_fusion.loc[
            key, f"Test {metric}"
        ] = f"{mean_score:.3f} (±{std_score:.3f})"

compare_fusion.to_csv(pf_table_fn, index=True, header=True)

## Perform statistical comparisons

# First obtain only the numbers to use for summary & statistical evaluation
test_performances = pd.DataFrame.from_dict(
    {
        k: v.iloc[
            list(v.index == "Average"), [x.startswith("Test ROC") for x in v.columns]
        ].values.ravel()
        for k, v in performances.items()
    }
)

for combo in combinations(test_performances.columns, 2):
    t, p = ttest_ind(test_performances[combo[0]], test_performances[combo[1]])
    print(f"Comparison - {combo[0]} vs. {combo[1]}: t={t:.2f}, p={p:.2E}")
