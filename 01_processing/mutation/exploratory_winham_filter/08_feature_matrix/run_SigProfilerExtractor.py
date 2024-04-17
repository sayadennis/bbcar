# pylint: disable=import-error

import sys

from SigProfilerExtractor import sigpro as sig

filter_type = sys.argv[1]

din = (
    "/projects/b1131/saya/bbcar/exploratory_winham_filter"
    f"/{filter_type}/07_predicted_somatic/vcfs"
)
dout = (
    "/projects/b1131/saya/bbcar/exploratory_winham_filter"
    f"/{filter_type}/08_feature_matrix/signature_results"
)

if __name__ == "__main__":
    sig.sigProfilerExtractor(
        input_type="vcf",
        output=dout,
        input_data=din,
        reference_genome="GRCh38",
        minimum_signatures=1,
        maximum_signatures=10,
        nmf_replicates=100,
        cpu=32,
    )
