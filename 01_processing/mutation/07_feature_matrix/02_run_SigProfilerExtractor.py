# pylint: disable=missing-module-docstring
# pylint: disable=import-error

import sys

from SigProfilerExtractor import sigpro as sig

din = sys.argv[1]
dout = sys.argv[2]
n_cpu = int(sys.argv[3])

print(f"Input directory: {din}")
print(f"Output directory: {dout}")
print(f"Number of CPUs: {n_cpu}")

if __name__ == "__main__":
    sig.sigProfilerExtractor(
        input_type="vcf",
        output=dout,
        input_data=din,
        reference_genome="GRCh38",
        minimum_signatures=1,
        maximum_signatures=12,
        cpu=n_cpu,
    )
