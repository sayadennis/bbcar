import glob

import numpy as np
import pandas as pd

########################################################################
#### Target from names of files in GATK segments that Gannon shared ####
########################################################################

din = "/projects/b1131/saya/bbcar/data/02b_cnv/01_gatk_analyzed_segments"
dout = "/projects/b1131/saya/bbcar/data/clinical"

sampleid_list = []
cc_list = []

for fn in glob.glob(din + "/*.csv"):
    sampleid = fn.split("/")[-1].split(".")[0].split("_")[0]
    cc = fn.split("/")[-1].split(".")[0].split("_")[1]
    sampleid_list.append(int(sampleid))
    if cc == "Tissue":
        cc_list.append(1)
    elif cc == "Control":
        cc_list.append(0)
    else:
        print(f"Unexpected value in case/control: {cc}")

label = pd.DataFrame(np.array(cc_list).T, index=sampleid_list, columns=["label"])
label.sort_index(axis=0, inplace=True)

label.to_csv(
    f"{dout}/bbcar_label_studyid_from_gatk_filenames.csv", header=True, index=True
)

#################################
#### Target from REDCap data ####
#################################

din = "/projects/b1131/saya/bbcar/data/00_raw"
dout = "/projects/b1131/saya/bbcar/data/clinical"

# read data
clin = pd.read_csv(f"{din}/BBCaRDatabaseNU09B2_DATA_2022-10-21_1440.csv")

# Study ID looks something like "BBC0001-500" -> change to integer IDs
clin["study_id_int"] = [int(x.split("-")[0][3:]) for x in clin.study_id.values]

# data includes lots of clinical data - only keep the case/control info for now
clin = clin[["study_id_int", "cc_id"]].set_index("study_id_int", drop=True)
clin.index.name = None

# REDCap mapping - obtained from https://redcap.nubic.northwestern.edu/redcap/redcap_v12.0.26/Design/data_dictionary_codebook.php?pid=3928 # pylint: disable=line-too-long
redcapdict = {1: "Case", 2: "Control"}

# Binarize target
clin.cc_id = clin.cc_id.map(redcapdict).map({"Case": 1, "Control": 0})

# Get the ID os patients whose sequencing data we have
record_ids = []
for source in ["germline", "tissue"]:
    record_ids += [int(x.split("/")[-1]) for x in glob.glob(f"{din}/{source}/[1-9]*")]

clin = clin.iloc[[i in record_ids for i in clin.index], :].astype(int)

# save
clin.to_csv(f"{dout}/bbcar_redcap_label_studyid.csv", header=True, index=True)

##############################
#### Summarize difference ####
##############################

# join the two targets by REDCap record ID
joined = label.merge(clin, how="inner", left_index=True, right_index=True)

# print the IDs where the labels don't match
disagree = []
for i in joined.index:
    if joined.loc[i, "label"] != joined.loc[i, "cc_id"]:
        disagree.append(i)  # non-matching between my previous & REDCap

print(joined.loc[disagree, :])
