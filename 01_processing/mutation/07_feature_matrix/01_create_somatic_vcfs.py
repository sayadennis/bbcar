import glob
import os

import pandas as pd

dn = "/projects/b1131/saya/new_bbcar/data/02a_mutation"

sompred = pd.read_csv(f"{dn}/07_predicted_somatic/nonmatched.csv")
sompred = list(sompred.iloc[sompred.somatic.values == 1, :].var_id.values)

tissue_vcfs = glob.glob(
    f"{dn}/02_variant_calls/tissue_only/*_DPfiltered_classicalAF.vcf"
)
tissue_sample_ids = [filename.split("/")[-1].split("_")[0] for filename in tissue_vcfs]

for sample_id in tissue_sample_ids:
    fin = f"{dn}/02_variant_calls/tissue_only/{sample_id}_DPfiltered_classicalAF.vcf"
    if not os.path.exists(f"{dn}/07_predicted_somatic/vcfs/"):
        os.makedirs(f"{dn}/07_predicted_somatic/vcfs/")
    #
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
