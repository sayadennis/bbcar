# pylint: disable=import-error
import sys

from SigProfilerExtractor import sigpro as sig

din = sys.argv[1]
dout = sys.argv[2]
n_cpu = int(sys.argv[3])

if __name__ == "__main__":
    sig.sigProfilerExtractor(
        input_type="seg:ASCAT",
        output=dout,
        input_data=f"{din}/all_ASCAT_segs_concat.txt",
        reference_genome="GRCh38",
        minimum_signatures=1,
        maximum_signatures=12,
        nmf_replicates=100,
        cpu=n_cpu,
    )
