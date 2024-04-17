import pandas as pd
from sklearn.impute import IterativeImputer

##########################################
#### Set input and output directories ####
##########################################

din = "/projects/b1131/saya/bbcar/data/02a_mutation/04_ml_features/03_freq_added"
dout = "/projects/b1131/saya/bbcar/data/02a_mutation/04_ml_features/04_imputed"

for pon_source in ["bbcar", "1000g"]:
    filename = f"features_annovar_bbcarfreq_{pon_source}pon.csv"

    ###################
    #### Load data ####
    ###################

    features = pd.read_csv(f"{din}/{filename}")

    ############################
    #### Numerical encoding ####
    ############################

    features.SIFT_pred = features.SIFT_pred.map({"D": 1, "T": 0})
    features.Polyphen2_HDIV_pred = features.Polyphen2_HDIV_pred.map(
        {"B": 0, "P": 1, "D": 2}
    )
    features.Polyphen2_HVAR_pred = features.Polyphen2_HVAR_pred.map(
        {"B": 0, "P": 1, "D": 2}
    )
    features.LRT_pred = features.LRT_pred.map({"N": 0, "U": 1, "D": 2})
    features.MutationTaster_pred = features.MutationTaster_pred.map(
        {"P": 0, "N": 1, "D": 2, "A": 3}
    )
    features.MutationAssessor_pred = features.MutationAssessor_pred.map(
        {"N": 0, "L": 1, "M": 2, "H": 3}
    )
    features.FATHMM_pred = features.FATHMM_pred.map({"T": 0, "D": 1})
    features.MetaSVM_pred = features.MetaSVM_pred.map({"T": 0, "D": 1})
    features.MetaLR_pred = features.MetaLR_pred.map({"T": 0, "D": 1})

    ####################
    #### Imputation ####
    ####################

    features_to_impute = [
        "AF",
        "ExAC_ALL",
        "SIFT_score",
        "SIFT_pred",
        "Polyphen2_HDIV_score",
        "Polyphen2_HDIV_pred",
        "Polyphen2_HVAR_score",
        "Polyphen2_HVAR_pred",
        "LRT_score",
        "LRT_pred",
        "MutationTaster_score",
        "MutationTaster_pred",
        "MutationAssessor_score",
        "MutationAssessor_pred",
        "FATHMM_score",
        "FATHMM_pred",
        "MetaSVM_score",
        "MetaSVM_pred",
        "MetaLR_score",
        "MetaLR_pred",
        "VEST3_score",
        "CADD_raw",
        "CADD_phred",
        "GERP++_RS",
        "phyloP20way_mammalian",
        "phyloP100way_vertebrate",
        "SiPhy_29way_logOdds",
    ]

    imp_median = IterativeImputer(
        random_state=0, initial_strategy="median", max_iter=20
    )
    imp_features = imp_median.fit_transform(features[features_to_impute])

    features[features_to_impute] = imp_features

    ## Save ##
    features.to_csv(
        f"{dout}/features_imputed_{pon_source}.csv", index=False, header=True
    )
