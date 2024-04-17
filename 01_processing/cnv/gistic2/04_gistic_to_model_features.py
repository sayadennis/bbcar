# pylint: disable=import-error

import os

import numpy as np
import pandas as pd
import process_functions

d01 = "/projects/b1131/saya/bbcar/data/02b_cnv/01_gatk_analyzed_segments"
d03 = "/projects/b1131/saya/bbcar/data/02b_cnv/03_gistic2_out_conf90"
d04 = "/projects/b1131/saya/bbcar/data/02b_cnv/04_cleaned_cnv"
dix = "/projects/b1131/saya/bbcar/data/02b_cnv/indices"

###############################
#### All samplies combined ####
###############################

cnv_dict = process_functions.generate_gistic_features(f"{d03}/all")

if not os.path.isdir(f"{d04}/all"):
    os.mkdir(f"{d04}/all")

cnv_dict["genecopy"].to_csv(
    f"{d04}/all/gene_copy_conf90_all.csv", header=True, index=True
)
cnv_dict["genethres"].to_csv(
    f"{d04}/all/gene_thres_conf90_all.csv", header=True, index=True
)
(cnv_dict["genecopy"] + 2).to_csv(
    f"{d04}/all/gene_copy_conf90_all_plus2.csv", header=True, index=True
)
(cnv_dict["genethres"] + 2).to_csv(
    f"{d04}/all/gene_thres_conf90_all_plus2.csv", header=True, index=True
)

cnv_dict["regcopy"].to_csv(
    f"{d04}/all/reg_copy_conf90_all.csv", header=True, index=True
)
cnv_dict["regthres"].to_csv(
    f"{d04}/all/reg_thres_conf90_all.csv", header=True, index=True
)
(cnv_dict["regcopy"] + 2).to_csv(
    f"{d04}/all/reg_copy_conf90_all_plus2.csv", header=True, index=True
)

cnv_dict["cytocopy"].to_csv(
    f"{d04}/all/cyto_copy_conf90_all.csv", header=True, index=True
)
cnv_dict["cytothres_amp"].to_csv(
    f"{d04}/all/cyto_thres_amp_conf90_all.csv", header=True, index=True
)
cnv_dict["cytothres_del"].to_csv(
    f"{d04}/all/cyto_thres_del_conf90_all.csv", header=True, index=True
)
(cnv_dict["cytocopy"] + 2).to_csv(
    f"{d04}/all/cyto_copy_conf90_all_plus2.csv", header=True, index=True
)

###############################################
#### All cases and all controls separately ####
###############################################

# create features for cases and controls
cnv_case = process_functions.generate_gistic_features(f"{d03}/all_cases")
cnv_cont = process_functions.generate_gistic_features(f"{d03}/all_controls")

# stack the gene-level features together
gene_copy = pd.concat((cnv_case["genecopy"], cnv_cont["genecopy"]), axis=0).sort_index(
    axis=0
)
gene_thres = pd.concat(
    (cnv_case["genethres"], cnv_cont["genethres"]), axis=0
).sort_index(axis=0)

# stack the cytoband-level features together
cyto_copy = (
    pd.concat((cnv_case["cytocopy"], cnv_cont["cytocopy"]), axis=0)
    .fillna(0.0)
    .sort_index(axis=0)
)
cyto_thres_amp = (
    pd.concat((cnv_case["cytothres_amp"], cnv_cont["cytothres_amp"]), axis=0)
    .fillna(0.0)
    .sort_index(axis=0)
)
cyto_thres_del = (
    pd.concat((cnv_case["cytothres_del"], cnv_cont["cytothres_del"]), axis=0)
    .fillna(0.0)
    .sort_index(axis=0)
)

if not os.path.isdir(f"{d04}/all_ccsep"):
    os.mkdir(f"{d04}/all_ccsep")

## save all the feature matrices
gene_copy.to_csv(f"{d04}/all_ccsep/gene_copy_conf90_ccsep.csv", header=True, index=True)
gene_thres.to_csv(
    f"{d04}/all_ccsep/gene_thres_conf90_ccsep.csv", header=True, index=True
)
(gene_thres + 2).to_csv(
    f"{d04}/all_ccsep/gene_thres_conf90_ccsep_plus2.csv", header=True, index=True
)

cnv_case["regcopy"].to_csv(
    f"{d04}/all_ccsep/reg_copy_conf90_cases.csv", header=True, index=True
)
cnv_cont["regcopy"].to_csv(
    f"{d04}/all_ccsep/reg_copy_conf90_controls.csv", header=True, index=True
)
cnv_case["regthres"].to_csv(
    f"{d04}/all_ccsep/reg_thres_conf90_cases.csv", header=True, index=True
)
cnv_cont["regthres"].to_csv(
    f"{d04}/all_ccsep/reg_thres_conf90_controls.csv", header=True, index=True
)

cyto_copy.to_csv(f"{d04}/all_ccsep/cyto_copy_conf90_ccsep.csv", header=True, index=True)
cyto_thres_amp.to_csv(
    f"{d04}/all_ccsep/cyto_thres_amp_conf90_ccsep.csv", header=True, index=True
)
cyto_thres_del.to_csv(
    f"{d04}/all_ccsep/cyto_thres_del_conf90_ccsep.csv", header=True, index=True
)

#########################################################
#### Training cases and training controls separately ####
#########################################################

cnv_case = process_functions.generate_gistic_features(f"{d03}/train_cases")
cnv_cont = process_functions.generate_gistic_features(f"{d03}/train_controls")

# stack the gene-level features together
gene_copy = pd.concat((cnv_case["genecopy"], cnv_cont["genecopy"]), axis=0).sort_index(
    axis=0
)
gene_thres = pd.concat(
    (cnv_case["genethres"], cnv_cont["genethres"]), axis=0
).sort_index(axis=0)

# stack the cytoband-level features together
cyto_copy = (
    pd.concat((cnv_case["cytocopy"], cnv_cont["cytocopy"]), axis=0)
    .fillna(0.0)
    .sort_index(axis=0)
)
cyto_thres_amp = (
    pd.concat((cnv_case["cytothres_amp"], cnv_cont["cytothres_amp"]), axis=0)
    .fillna(0.0)
    .sort_index(axis=0)
)
cyto_thres_del = (
    pd.concat((cnv_case["cytothres_del"], cnv_cont["cytothres_del"]), axis=0)
    .fillna(0.0)
    .sort_index(axis=0)
)

if not os.path.isdir(f"{d04}/train_ccsep"):
    os.mkdir(f"{d04}/train_ccsep")

## save all the feature matrices
gene_copy.to_csv(
    f"{d04}/train_ccsep/gene_copy_conf90_ccsep.csv", header=True, index=True
)
gene_thres.to_csv(
    f"{d04}/train_ccsep/gene_thres_conf90_ccsep.csv", header=True, index=True
)
(gene_thres + 2).to_csv(
    f"{d04}/train_ccsep/gene_thres_conf90_ccsep_plus2.csv", header=True, index=True
)

cnv_case["regcopy"].to_csv(
    f"{d04}/train_ccsep/reg_copy_conf90_cases.csv", header=True, index=True
)
cnv_cont["regcopy"].to_csv(
    f"{d04}/train_ccsep/reg_copy_conf90_controls.csv", header=True, index=True
)
cnv_case["regthres"].to_csv(
    f"{d04}/train_ccsep/reg_thres_conf90_cases.csv", header=True, index=True
)
cnv_cont["regthres"].to_csv(
    f"{d04}/train_ccsep/reg_thres_conf90_controls.csv", header=True, index=True
)

cyto_copy.to_csv(
    f"{d04}/train_ccsep/cyto_copy_conf90_ccsep.csv", header=True, index=True
)
cyto_thres_amp.to_csv(
    f"{d04}/train_ccsep/cyto_thres_amp_conf90_ccsep.csv", header=True, index=True
)
cyto_thres_del.to_csv(
    f"{d04}/train_ccsep/cyto_thres_del_conf90_ccsep.csv", header=True, index=True
)

###################################################################################
#### (Training cases + test all) and (training controls + test all) separately ####
###################################################################################

cnv_case = process_functions.generate_gistic_features(f"{d03}/train_cases_test_all")
cnv_cont = process_functions.generate_gistic_features(f"{d03}/train_controls_test_all")

# stack the gene-level features together
gene_copy = pd.concat((cnv_case["genecopy"], cnv_cont["genecopy"]), axis=0).sort_index(
    axis=0
)
gene_thres = pd.concat(
    (cnv_case["genethres"], cnv_cont["genethres"]), axis=0
).sort_index(axis=0)

# stack the cytoband-level features together
cyto_copy = (
    pd.concat((cnv_case["cytocopy"], cnv_cont["cytocopy"]), axis=0)
    .fillna(0.0)
    .sort_index(axis=0)
)
cyto_thres_amp = (
    pd.concat((cnv_case["cytothres_amp"], cnv_cont["cytothres_amp"]), axis=0)
    .fillna(0.0)
    .sort_index(axis=0)
)
cyto_thres_del = (
    pd.concat((cnv_case["cytothres_del"], cnv_cont["cytothres_del"]), axis=0)
    .fillna(0.0)
    .sort_index(axis=0)
)

if not os.path.isdir(f"{d04}/train_ccsep_test_all"):
    os.mkdir(f"{d04}/train_ccsep_test_all")

## save all the feature matrices
gene_copy.to_csv(
    f"{d04}/train_ccsep_test_all/gene_copy_conf90_ccsep.csv", header=True, index=True
)
gene_thres.to_csv(
    f"{d04}/train_ccsep_test_all/gene_thres_conf90_ccsep.csv", header=True, index=True
)
(gene_thres + 2).to_csv(
    f"{d04}/train_ccsep_test_all/gene_thres_conf90_ccsep_plus2.csv",
    header=True,
    index=True,
)

cnv_case["regcopy"].to_csv(
    f"{d04}/train_ccsep_test_all/reg_copy_conf90_cases.csv", header=True, index=True
)
cnv_cont["regcopy"].to_csv(
    f"{d04}/train_ccsep_test_all/reg_copy_conf90_controls.csv", header=True, index=True
)
cnv_case["regthres"].to_csv(
    f"{d04}/train_ccsep_test_all/reg_thres_conf90_cases.csv", header=True, index=True
)
cnv_cont["regthres"].to_csv(
    f"{d04}/train_ccsep_test_all/reg_thres_conf90_controls.csv", header=True, index=True
)

cyto_copy.to_csv(
    f"{d04}/train_ccsep_test_all/cyto_copy_conf90_ccsep.csv", header=True, index=True
)
cyto_thres_amp.to_csv(
    f"{d04}/train_ccsep_test_all/cyto_thres_amp_conf90_ccsep.csv",
    header=True,
    index=True,
)
cyto_thres_del.to_csv(
    f"{d04}/train_ccsep_test_all/cyto_thres_del_conf90_ccsep.csv",
    header=True,
    index=True,
)


##########################################################################
##########################################################################
######## Transform the test sets to generate #############################
######## final model features with both training and testing data ########
##########################################################################
##########################################################################

## Get test-set index
all_ix = pd.read_csv(f"{dix}/index_dictionary.csv", index_col=0)
all_studyid = all_ix["studyid"]

test_ix = pd.read_csv(f"{dix}/test_ix.csv", header=None)[0]  # this is index 0-200
test_studyid = list(
    pd.read_csv(f"{dix}/test_studyid.csv", header=None)[0]
)  # this is the study ID that matches the GATK file names

## Get cytoband coordinates and re-format to use it in "transform_gatk_to_reg_matrix"
cyto_coords = pd.read_csv(
    "/projects/b1122/saya/ref_cytobands/cytoband_hg38.tsv", sep="\t"
)
cyto_coords = cyto_coords.iloc[~pd.isnull(cyto_coords["name"].values), :].reset_index(
    drop=True
)  # only select cytobands with names (characterized)
cyto_regions = pd.DataFrame(
    {
        "chrom": cyto_coords["#chrom"],
        "start": cyto_coords["chromStart"],
        "end": cyto_coords["chromEnd"],
    }
)
cyto_regions.index = [
    cyto_coords["#chrom"].iloc[i].split("chr")[1] + cyto_coords["name"].iloc[i]
    for i in cyto_coords.index
]

## Get gene coordinates and re-format to use it in "transform_gatk_to_reg_matrix"
gene_coords = pd.read_csv(
    "/projects/b1122/saya/ref_gene_annotations/GRCh38_gencode.v27.refFlat.txt",
    sep="\t",
    header=None,
)
chroms = list(cyto_regions["chrom"].unique())
gene_coords = gene_coords.iloc[
    [x in chroms for x in gene_coords[2].values], :
].reset_index(
    drop=True
)  # only select normal chromosomes
gene_coords = gene_coords[~gene_coords[0].duplicated(keep="first")].reset_index(
    drop=True
)  # only keep one instance of each gene
gene_regions = pd.DataFrame(
    {"chrom": gene_coords[2], "start": gene_coords[4], "end": gene_coords[5]}
)
gene_regions.index = list(gene_coords[0].values)

###############################################
#### All cases and all controls separately ####
###############################################

case_regions = process_functions.lesions_to_structured_regions(
    f"{d03}/all_cases/all_lesions.conf_90.txt"
)
control_regions = process_functions.lesions_to_structured_regions(
    f"{d03}/all_controls/all_lesions.conf_90.txt"
)

ol_mx = process_functions.generate_ovelap_mx(case_regions, control_regions)

ol_thres = 1e2  # minimum amount of overlap between regions to consider "overlapping"
case_unique_regions = list(ol_mx.index[ol_mx.sum(axis=1).values < ol_thres])
control_unique_regions = list(ol_mx.columns[ol_mx.sum(axis=0).values < ol_thres])

# create a reg_table with these unique regions
reg_table = pd.concat(
    (
        case_regions.iloc[[x in case_unique_regions for x in case_regions.index]],
        control_regions.iloc[
            [x in control_unique_regions for x in control_regions.index]
        ],
    ),
    axis=0,
)
# reindex with case/control info to make indices unique
reg_table.index = ["Case " + x for x in case_unique_regions] + [
    "Control " + x for x in control_unique_regions
]

allsample_mx = process_functions.transform_gatk_to_reg_matrix(
    reg_table=reg_table, gatk_dir=d01
)
allsample_mx.to_csv(f"{d04}/all_ccsep/reg_copy_unique.csv", header=True, index=True)


#########################################################
#### Training cases and training controls separately ####
#########################################################

## Cytoband-level features - record the amp/del with largest absolute value within each cytoband
train_cytocopy = pd.read_csv(
    f"{d04}/train_ccsep/cyto_copy_conf90_ccsep.csv", index_col=0
)  # read training matrix
aber_cyto_regions = cyto_regions.iloc[
    [x in list(train_cytocopy.columns) for x in cyto_regions.index], :
]

test_cytocopy = process_functions.transform_gatk_to_reg_matrix(
    aber_cyto_regions, d01, sampleset=test_studyid
)

ccsep_cytocopy = pd.concat((train_cytocopy, test_cytocopy), axis=0).sort_index(axis=0)
ccsep_cytocopy.to_csv(f"{d04}/train_ccsep/cyto_copy.csv", header=True, index=True)

## Gene-level features - if sample has amp/del seg overlapping gene, enter the log copy number
# Get gene feature list from cleaned training set
train_genecopy = pd.read_csv(
    f"{d04}/train_ccsep/gene_copy_conf90_ccsep.csv", index_col=0
)  # read training matrix
gene_regions = gene_regions.iloc[
    [x in list(train_genecopy.columns) for x in gene_regions.index], :
]

test_genecopy = process_functions.transform_gatk_to_reg_matrix(
    gene_regions, d01, sampleset=list(test_studyid)
)

ccsep_genecopy = pd.concat((train_genecopy, test_genecopy), axis=0).sort_index(axis=0)
ccsep_genecopy.to_csv(f"{d04}/train_ccsep/gene_copy.csv", header=True, index=True)

## Region-level features - get case-/control-unique regions, record CN value of overlapping seg
# read case and control regions
case_regions = process_functions.lesions_to_structured_regions(
    f"{d03}/train_cases/all_lesions.conf_90.txt"
)
control_regions = process_functions.lesions_to_structured_regions(
    f"{d03}/train_controls/all_lesions.conf_90.txt"
)
# compute overlaps
ol_mx = process_functions.generate_ovelap_mx(case_regions, control_regions)
# find number of overlapping regions for case & control respectively (unnecessary)
olsize_control = np.sum(ol_mx.sum(axis=0).values > ol_thres)
olsize_case = np.sum(ol_mx.sum(axis=1).values > ol_thres)
# get a list of case-unique and control-unique region names
case_unique_regions = list(ol_mx.index[ol_mx.sum(axis=1).values < ol_thres])
control_unique_regions = list(ol_mx.columns[ol_mx.sum(axis=0).values < ol_thres])
# create a reg_table with these unique regions
reg_table = pd.concat(
    (
        case_regions.iloc[[x in case_unique_regions for x in case_regions.index]],
        control_regions.iloc[
            [x in control_unique_regions for x in control_regions.index]
        ],
    ),
    axis=0,
)
# reindex with case/control info to make row-indices unique
reg_table.index = ["Case " + x for x in case_unique_regions] + [
    "Control " + x for x in control_unique_regions
]
# save
allsample_mx = process_functions.transform_gatk_to_reg_matrix(
    reg_table=reg_table, gatk_dir=d01
)
allsample_mx.to_csv(f"{d04}/train_ccsep/reg_copy_unique.csv", header=True, index=True)


###################################################################################
#### (Training cases + test all) and (training controls + test all) separately ####
###################################################################################

## Cytoband-level features
cytocopy = pd.read_csv(
    f"{d04}/train_ccsep_test_all/cyto_copy_conf90_ccsep.csv", index_col=0
)
cytocopy_transformed = cytocopy.iloc[[x not in test_studyid for x in cytocopy.index], :]
cytocopy_test = cytocopy.iloc[[x in test_studyid for x in cytocopy.index], :]
for studyid in cytocopy_test.index.unique():
    subset = cytocopy_test.iloc[cytocopy_test.index == studyid, :]
    cytocopy_transformed.loc[studyid, :] = list(subset.mean())

cytocopy_transformed.sort_index(axis=0).to_csv(
    f"{d04}/train_ccsep_test_all/cyto_copy.csv", header=True, index=True
)

## Gene-level features
genecopy = pd.read_csv(
    f"{d04}/train_ccsep_test_all/gene_copy_conf90_ccsep.csv", index_col=0
)
genecopy_transformed = genecopy.iloc[[x not in test_studyid for x in genecopy.index], :]
genecopy_test = genecopy.iloc[[x in test_studyid for x in genecopy.index], :]
for studyid in genecopy_test.index.unique():
    subset = genecopy_test.iloc[genecopy_test.index == studyid, :]
    genecopy_transformed.loc[studyid, :] = list(subset.mean())

genecopy_transformed.sort_index(axis=0).to_csv(
    f"{d04}/train_ccsep_test_all/gene_copy.csv", header=True, index=True
)

## Region-level features
# read case and control regions
case_regions = process_functions.lesions_to_structured_regions(
    f"{d03}/train_cases_test_all/all_lesions.conf_90.txt"
)
control_regions = process_functions.lesions_to_structured_regions(
    f"{d03}/train_controls_test_all/all_lesions.conf_90.txt"
)
# compute overlaps
ol_mx = process_functions.generate_ovelap_mx(case_regions, control_regions)
# find number of overlapping regions for case & control respectively (unnecessary)
olsize_control = np.sum(ol_mx.sum(axis=0).values != 0)
olsize_case = np.sum(ol_mx.sum(axis=1).values != 0)
# get a list of case-unique and control-unique region names
case_unique_regions = list(ol_mx.index[ol_mx.sum(axis=1).values == 0])
control_unique_regions = list(ol_mx.columns[ol_mx.sum(axis=0).values == 0])
# create a reg_table with these unique regions
reg_table = pd.concat(
    (
        case_regions.iloc[[x in case_unique_regions for x in case_regions.index]],
        control_regions.iloc[
            [x in control_unique_regions for x in control_regions.index]
        ],
    ),
    axis=0,
)
# reindex with case/control info to make indices unique
reg_table.index = ["Case " + x for x in case_unique_regions] + [
    "Control " + x for x in control_unique_regions
]
# save
allsample_mx = process_functions.transform_gatk_to_reg_matrix(
    reg_table=reg_table, gatk_dir=d01
)
allsample_mx.to_csv(
    f"{d04}/train_ccsep_test_all/reg_copy_unique.csv", header=True, index=True
)
