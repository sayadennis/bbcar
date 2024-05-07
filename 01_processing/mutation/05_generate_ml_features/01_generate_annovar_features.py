# pylint: disable=duplicate-code

import os
import sys

import numpy as np
import pandas as pd

###########################################
#### Get input arguments (sample info) ####
###########################################

patient_id = sys.argv[1]
print(patient_id)
source = sys.argv[2]
print(source)

##########################################
#### Set input and output directories ####
##########################################

din = f"/projects/b1131/saya/new_bbcar/data/02a_mutation/03_annotated_variants/annovar/{source}"
dout = "/projects/b1131/saya/new_bbcar/data/02a_mutation/04_ml_features/individual_annovar_features"

if not os.path.exists(dout):
    os.makedirs(dout, exist_ok=True)

#######################################
#### Set variable names to collect ####
#######################################

with open("av_features.txt", "r", encoding="utf-8") as f:
    av_features = [x.strip() for x in f.readlines()]

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
] + av_features

##########################
#### Collect features ####
##########################

## Define feature matrix where I can record/append data
features = pd.DataFrame(columns=variables)

## Iterate through lines of VCF
vcf = f"{din}/{patient_id}.hg38_multianno.vcf"

with open(vcf, "r") as f:
    lines = f.readlines()

# loop through lines
for line in lines:
    if not line.startswith("#"):
        var_feature = pd.DataFrame(index=[0], columns=variables)
        # define elements of variants and their annotations
        var_elements = line.split()
        chrom = var_elements[0]
        pos = var_elements[1]
        ref = var_elements[3]
        alt = var_elements[4]
        info = var_elements[7]
        # create variant ID
        var_feature.loc[0, "var_id"] = f"{chrom}_{pos}_{ref}_{alt}"
        var_feature.loc[0, "source"] = source
        var_feature.loc[0, "sample_id"] = patient_id
        # get necessary information
        info = info.split(";")
        for variable in variables:
            if variable == "AF":
                af_locs = np.where([anno_item.startswith("AF=") for anno_item in info])
                afs = [x.split("=")[1] for x in np.array(info)[af_locs]]
                afs = [item for item in afs if item != "."]
                afs = [
                    item.split(",")[0] if item == "0.5,0.5" else item for item in afs
                ]
                if len(afs) > 0:
                    af = np.mean([float(item) for item in afs])
                    var_feature.loc[0, variable] = af
            elif variable in ["var_id", "source", "sample_id"]:
                continue
            else:
                feature_value = [
                    item
                    for (item, filterval) in zip(
                        info,
                        [anno_item.startswith(f"{variable}=") for anno_item in info],
                    )
                    if filterval
                ]
                feature_value = feature_value[0].split("=")[-1]
                if feature_value.replace(".", "").replace("-", "").isnumeric():
                    var_feature.loc[0, variable] = float(feature_value)
                else:
                    var_feature.loc[0, variable] = (
                        feature_value if feature_value != "." else None
                    )
        # record result in feature matrix
        features = pd.concat((features, var_feature), ignore_index=True)

features.to_csv(f"{dout}/{source}_{patient_id}_annovar_features.csv", index=False)
