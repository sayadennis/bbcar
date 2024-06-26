import os

import pandas as pd

din = "/projects/b1131/saya/new_bbcar/data/02b_cnv/signatures/03_ASCAT_obj"
dout = din

tissue_sampleid_fn = (
    "/projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_tissue.txt"
)
germline_sampleid_fn = (
    "/projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_germline.txt"
)

with open(tissue_sampleid_fn, "r", encoding="utf-8") as f:
    tissue_sampleids = set(int(x.strip()) for x in f.readlines())

with open(germline_sampleid_fn, "r", encoding="utf-8") as f:
    germline_sampleids = set(int(x.strip()) for x in f.readlines())

matched_sampleids = list(germline_sampleids.intersection(tissue_sampleids))

segs = pd.DataFrame()

for sampleid in matched_sampleids:
    if os.path.exists(f"{din}/{sampleid}_tissue.segments.txt"):
        sample_segs = pd.read_csv(f"{din}/{sampleid}_tissue.segments.txt", sep="\t")
        segs = pd.concat((segs, sample_segs), axis=0)
    else:
        print(f"Tissue segments does not exist for sample {sampleid}")

segs["chr"] = [x[3:] for x in segs["chr"]]
segs["sample"] = [x.split("_")[0] for x in segs["sample"]]

segs.to_csv(f"{dout}/all_ASCAT_segs_concat.txt", sep="\t", index=False, header=True)
