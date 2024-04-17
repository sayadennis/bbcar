# pylint: disable=duplicate-code

import pandas as pd

##########################################
#### Set input and output directories ####
##########################################

din = "/projects/b1131/saya/bbcar/data/02a_mutation/04_ml_features/02_concat_annovar_features"
dout = "/projects/b1131/saya/bbcar/data/02a_mutation/04_ml_features/03_freq_added"

############################################
#### Set empty dataframe to concatenate ####
############################################

with open("av_features.txt", "r", encoding="utf-8") as f:
    av_features = [x.strip() for x in f.readlines()]

variables = ["var_id", "source", "sample_id"] + av_features

for pon_source in ["1000g", "bbcar"]:
    # load the concatenated annovar features
    features = pd.read_csv(f"{din}/annovar_features_all_{pon_source}pon.csv")
    # create binary matrix where index is var_id and column is sample_id
    matrix = features[["var_id", "sample_id"]].pivot_table(
        index="var_id", columns="sample_id", aggfunc=lambda x: 1, fill_value=0
    )
    # calculate the frequency among bbcar samples and save it in a vector
    bbcar_freq = matrix.sum(axis=1) / matrix.shape[1]
    # add this vector as a column to the existing features matrix
    features["bbcar_freq"] = bbcar_freq.loc[features.var_id.values].values
    # save the matrix
    features.to_csv(
        f"{dout}/features_annovar_bbcarfreq_{pon_source}pon.csv", index=False
    )
