import os

import pandas as pd

dn = "/projects/b1131/saya/new_bbcar/data/02a_mutation"

sompred = pd.read_csv(f"{dn}/07_predicted_somatic/nonmatched.csv")
sompred = list(sompred.iloc[sompred.somatic.values == 1, :].var_id.values)

with open(
    "/projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_tissue.txt", "r"
) as f:
    tissue_ids = [x.strip() for x in f.readlines()]

with open(
    "/projects/b1131/saya/new_bbcar/jobarray_args/patient_ids_germline.txt", "r"
) as f:
    germline_ids = [x.strip() for x in f.readlines()]

if not os.path.exists(f"{dn}/07_predicted_somatic/vcfs/"):
    os.makedirs(f"{dn}/07_predicted_somatic/vcfs/")

for sample_id in list(set(tissue_ids) - set(germline_ids)):
    fin = f"{dn}/02_variant_calls/tissue_only/{sample_id}_DPfiltered_classicalAF.vcf"
    fout = f"{dn}/07_predicted_somatic/vcfs/{sample_id}_somatic.vcf"
    # read the original calls
    with open(fin, "r") as f:
        lines = f.readlines()
    # write the filtered calls
    with open(fout, "w") as f:
        for line in lines:
            if line.startswith("#"):
                f.write(line)
            else:
                chrom = line.split()[0]
                pos = line.split()[1]
                ref = line.split()[3]
                alt = line.split()[4]
                var_id = f"{chrom}_{pos}_{ref}_{alt}"
                if var_id in sompred:
                    f.write(line)
