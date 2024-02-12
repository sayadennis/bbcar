# pylint: disable=missing-module-docstring

from itertools import combinations

import pandas as pd
from scipy.stats import ttest_ind

din = "/projects/b1131/saya/bbcar/model_interpretations/breast_cancer_prediction"

performances = {
    # "no fusion mutation" : pd.read_csv(f"{din}/no_fusion_mutation_orig.csv", index_col=0),
    "no fusion cnv": pd.read_csv(f"{din}/no_fusion_cnv_orig.csv", index_col=0),
    "early fusion": pd.read_csv(f"{din}/early_fusion_orig.csv", index_col=0),
    "mid fusion supervised": pd.read_csv(
        f"{din}/mid_fusion_supervised.csv", index_col=0
    ),
    "mid fusion unsupervised": pd.read_csv(
        f"{din}/mid_fusion_unsupervised.csv", index_col=0
    ),
}

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
