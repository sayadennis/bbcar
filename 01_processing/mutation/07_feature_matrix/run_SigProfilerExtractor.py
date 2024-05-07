# pylint: disable=missing-module-docstring
# pylint: disable=import-error

import pandas as pd
from SigProfilerExtractor import sigpro as sig

din = "/projects/b1131/saya/bbcar/data/02a_mutation/07_predicted_somatic/vcfs_matched"
dout = (
    "/projects/b1131/saya/bbcar/data/02a_mutation"
    "/08_feature_matrix/signature_results_matched_gpu"
)

if __name__ == "__main__":
    sig.sigProfilerExtractor(
        input_type="vcf",
        output=dout,
        input_data=din,
        reference_genome="GRCh38",
        exome=True,
        minimum_signatures=1,
        maximum_signatures=15,
        nmf_replicates=100,
        cpu=32,
        nmf_init="nndsvd_min",
        seeds="/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix/SigPro_Seeds.txt",
        gpu=True,
    )

data = pd.read_csv(f"{dout}/SBS96/Samples.txt", sep="\t", index_col=0).T
data.index = [x.split("_")[0] for x in data.index]

data.to_csv(f"{dout}/sbs_96_original_per_sample.csv")
