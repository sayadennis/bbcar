import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

################################
#### Get signature exposure ####
################################

sig_dir = "/projects/b1131/saya/new_bbcar/data/02b_cnv/signatures/04_signatures"
fin = f"{sig_dir}/CNV48/Samples.txt"

cn_features = pd.read_csv(fin, sep="\t")

cosmic_features = pd.read_csv(
    (
        f"{sig_dir}/CNV48/Suggested_Solution/COSMIC_CNV48_Decomposed_Solution"
        f"/Signatures/COSMIC_CNV48_Signatures.txt"
    ),
    sep="\t",
)

exposures = pd.DataFrame(
    np.dot(cn_features.iloc[:, 1:].T, cosmic_features.iloc[:, 1:]),
    index=cn_features.columns[1:],
    columns=cosmic_features.columns[1:],
)

######################
#### Get metadata ####
######################

meta = pd.read_csv("/projects/b1131/saya/new_bbcar/meta.csv")
case_ids = meta.iloc[meta.label.values == 1, :].patient_id.unique()

#################################
#### Compare exposure levels ####
#################################

for sig in exposures.columns:
    case_exposures = exposures.iloc[[x in case_ids for x in exposures.index], :][
        sig
    ].values.ravel()
    control_exposures = exposures.iloc[[x not in case_ids for x in exposures.index], :][
        sig
    ].values.ravel()
    t, p = ttest_ind(case_exposures, control_exposures)
    print(f"Signature {sig}: t={t:.2f}, p={p:.4f}")
