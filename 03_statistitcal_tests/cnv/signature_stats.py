import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

################################
#### Get signature exposure ####
################################

sig_dir = "/projects/b1131/saya/new_bbcar/data/02b_cnv/signatures/04_signatures"
fin = f"{sig_dir}/CNV48/Samples.txt"

cn_features = pd.read_csv(fin, sep="\t", index_col=0).T
cn_features.index = cn_features.index.astype(int)

denovo_sig = pd.read_csv(
    (
        f"{sig_dir}/CNV48/Suggested_Solution/CNV48_De-Novo_Solution"
        f"/Signatures/CNV48_De-Novo_Signatures.txt"
    ),
    sep="\t",
    index_col=0,
)

cosmic_sig = pd.read_csv(
    (
        f"{sig_dir}/CNV48/Suggested_Solution/COSMIC_CNV48_Decomposed_Solution"
        f"/Signatures/COSMIC_CNV48_Signatures.txt"
    ),
    sep="\t",
    index_col=0,
)

cosmic_exposures = pd.DataFrame(
    np.dot(cn_features, cosmic_sig),
    index=cn_features.index,
    columns=cosmic_sig.columns,
)

denovo_exposures = pd.DataFrame(
    np.dot(cn_features, denovo_sig),
    index=cn_features.index,
    columns=denovo_sig.columns,
)

######################
#### Get metadata ####
######################

meta = pd.read_csv("/projects/b1131/saya/new_bbcar/meta.csv")
case_ids = meta.iloc[meta.label.values == 1, :].patient_id.unique()

#################################
#### Compare exposure levels ####
#################################

exposures = {
    "denovo": denovo_exposures,
    "cosmic": cosmic_exposures,
}

for sigtype, df in exposures.items():
    print(f"#### {sigtype.upper()} ####")
    for sig in df.columns:
        case_exposures = df.iloc[[x in case_ids for x in df.index], :][
            sig
        ].values.ravel()
        control_exposures = df.iloc[[x not in case_ids for x in df.index], :][
            sig
        ].values.ravel()
        t, p = ttest_ind(case_exposures, control_exposures)
        print(f"Signature {sig}: t={t:.2f}, p={p:.4f}")
