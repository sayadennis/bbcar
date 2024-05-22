import numpy as np
import pandas as pd

din = "/projects/b1131/saya/new_bbcar/data/02a_mutation/08_feature_matrix/signature_results"
dout = "/projects/b1131/saya/new_bbcar/data/02a_mutation/08_feature_matrix"

for sigtype in ["SBS96", "DBS78", "ID83"]:
    sample_features = pd.read_csv(
        f"{din}/{sigtype}/Samples.txt", sep="\t", index_col=0
    ).T
    sample_features.index = [int(x.split("_")[0]) for x in sample_features.index]
    sample_features.to_csv(f"{dout}/raw_{sigtype}_features.csv")
    opt_num_sig = pd.read_csv(
        (
            f"{din}/{sigtype}/Suggested_Solution/{sigtype}_De-Novo_Solution/"
            f"Signatures/{sigtype}_De-Novo_Signatures.txt"
        ),
        sep="\t",
        index_col=0,
    ).shape[1]
    print(f"Using {opt_num_sig} signatures for {sigtype}")
    #### De Novo Signatures ####
    opt_solution = pd.read_csv(
        (
            f"{din}/{sigtype}/All_Solutions/{sigtype}_{opt_num_sig}_Signatures/"
            f"Signatures/{sigtype}_S{opt_num_sig}_Signatures.txt"
        ),
        sep="\t",
        index_col=0,
    )
    denovo_sigs_per_sample = np.dot(
        sample_features, opt_solution
    )  # realized that I had duplicate VCFs for matched samples - worknig on this now
    pd.DataFrame(
        denovo_sigs_per_sample,
        index=sample_features.index,
        columns=opt_solution.columns,
    ).to_csv(f"{dout}/denovo_signature_per_sample_{sigtype}.csv")
    #### Decomposed to COSMIC Signatures ####
    decomposed_sig = pd.read_csv(
        (
            f"{din}/{sigtype}/Suggested_Solution/COSMIC_{sigtype}_Decomposed_Solution"
            f"/Signatures/COSMIC_{sigtype}_Signatures.txt"
        ),
        sep="\t",
        index_col=0,
    )
    decomposed_sigs_per_sample = pd.DataFrame(
        np.dot(sample_features, decomposed_sig),
        index=sample_features.index,
        columns=decomposed_sig.columns,
    )
    decomposed_sigs_per_sample.to_csv(
        f"{dout}/cosmic_signature_per_sample_{sigtype}.csv"
    )
