# pylint: disable=unused-import

import glob
import os

import pandas as pd
from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer
from sklearn.model_selection import train_test_split

##########################################
#### Set input and output directories ####
##########################################

din = "/projects/b1131/saya/new_bbcar/data/02a_mutation/04_ml_features/individual_annovar_features"
dout = "/projects/b1131/saya/new_bbcar/data/02a_mutation/04_ml_features"
dix = "/projects/b1131/saya/new_bbcar/data/02a_mutation/04_ml_features/somatic_pred_ix"

if not os.path.exists(dix):
    os.makedirs(dix, exist_ok=True)

############################################
#### Set empty dataframe to concatenate ####
############################################

variables = [
    "var_id",
    "source",
    "sample_id",
    "Func.refGene",
    "Gene.refGene",
    "ExonicFunc.refGene",
    "Func.knownGene",
    "Gene.knownGene",
    "GeneDetail.knownGene",
    "ExonicFunc.knownGene",
    "Func.ensGene",
    "Gene.ensGene",
    "GeneDetail.ensGene",
    "ExonicFunc.ensGene",
    "AF",
    "avsnp150",
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

features = pd.DataFrame(columns=variables)

for fn in glob.glob(f"{din}/*_annovar_features.csv"):
    single_sample = pd.read_csv(fn)
    features = pd.concat((features, single_sample), ignore_index=True)

# # BELOW LINE MIGHT NOT BE NECESSARY if the weird '0.5,0.5' doesn't appear in AF column
# features.AF = features.AF.map({'0.5,0.5':0.5}).astype(float)

#### Change contents of the avsnp150 column to be binary ####
features["avsnp150"].fillna(0, inplace=True)
features["avsnp150"].iloc[features.avsnp150.values != 0] = 1

###########################################################
#### Calculate mutation frequency within BBCaR smaples ####
###########################################################

# create binary matrix where index is var_id and column is sample_id
matrix = features[["var_id", "sample_id"]].pivot_table(
    index="var_id", columns="sample_id", aggfunc=lambda x: 1, fill_value=0
)
# calculate the frequency among bbcar samples and save it in a vector
bbcar_freq = matrix.sum(axis=1) / matrix.shape[1]
# add this vector as a column to the existing features matrix
features["bbcar_freq"] = bbcar_freq.loc[features.var_id.values].values

###############################
#### Impute missing values ####
###############################

## Numerical encoding

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

## Imputation

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

imp_median = IterativeImputer(random_state=0, initial_strategy="median", max_iter=20)
imp_features = imp_median.fit_transform(features[features_to_impute])

features[features_to_impute] = imp_features

## Save ##
features.to_csv(f"{dout}/features_imputed.csv", index=False, header=True)

##########################
#### Train test split ####
##########################

data = features
data["somatic"] = (data.source == "tissue_normal").astype(int)

with open(f"{dout}/tissue_normal_var_id.txt", "w") as f:
    for var_id in list(data.iloc[data.source.values == "tissue_normal", :].var_id):
        f.write(f"{var_id}\n")

with open(f"{dout}/tissue_only_var_id.txt", "w") as f:
    for var_id in list(data.iloc[data.source.values == "tissue_only", :].var_id):
        f.write(f"{var_id}\n")

with open(f"{dout}/germline_var_id.txt", "w") as f:
    for var_id in list(data.iloc[data.source.values == "germline_only", :].var_id):
        f.write(f"{var_id}\n")

variables = [
    "var_id",
    # 'AF',
    "avsnp150",
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
    "somatic",
]

Xy_nonmatched = (
    data.iloc[data.source.values == "tissue_only", :][variables]
    .drop_duplicates(ignore_index=True)
    .set_index("var_id", drop=True)
)
X_nonmatched = Xy_nonmatched.iloc[:, :-1]

Xy_matched = (
    data.iloc[data.source.values != "tissue_only", :][variables]
    .drop_duplicates(ignore_index=True)
    .set_index("var_id", drop=True)
)
X_matched = Xy_matched.iloc[:, :-1]
y_matched = Xy_matched.iloc[:, -1]

# Subsample to decrease model training time
X_matched = X_matched.sample(n=int(2e5), replace=False, random_state=9)
y_matched = y_matched.loc[X_matched.index, :]

## Create train and test indices

train_ix, test_ix = train_test_split(
    X_matched.index, test_size=0.2, random_state=43, shuffle=True
)

## Save

X_matched.to_csv(f"{dout}/input_matched.csv", index=True, header=True)
y_matched.to_csv(f"{dout}/target_matched.csv", index=True, header=True)

X_nonmatched.to_csv(f"{dout}/input_nonmatched.csv", index=True, header=True)

pd.DataFrame(index=train_ix).to_csv(f"{dix}/train_index.txt", header=False)
pd.DataFrame(index=test_ix).to_csv(f"{dix}/test_index.txt", header=False)
