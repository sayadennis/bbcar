import os
import numpy as np
import pandas as pd
import glob

sampleid_list = []
cc_list = []

for fn in glob.glob("/projects/b1122/saya/bbcar_project/gatk_analyzed_segments/*.csv"):
    sampleid = fn.split("/")[-1].split(".")[0] # .split("_")[0]
    cc = fn.split("/")[-1].split(".")[0].split("_")[1]
    sampleid_list.append(sampleid)
    if cc=="Tissue":
        cc_list.append(1)
    elif cc=="Control":
        cc_list.append(0)
    else:
        print("Unexpected value in case/control: {}".format(cc))

label = pd.DataFrame(np.array(cc_list).T, index=sampleid_list, columns=["label"])

label.to_csv("/projects/b1122/saya/bbcar_project/cnv_matrix/bbcar_label.csv", header=True, index=True)
