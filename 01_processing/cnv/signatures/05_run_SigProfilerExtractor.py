# pylint: disable=import-error

from SigProfilerExtractor import sigpro as sig

din = "/projects/b1131/saya/bbcar/data/02b_cnv/signatures/03_ASCAT_obj/tissue_normal"
dout = "/projects/b1131/saya/bbcar/data/02b_cnv/signatures/04_signatures"

if __name__ == "__main__":
    sig.sigProfilerExtractor(
        input_type="seg:ASCAT",
        output=dout,
        input_data=f"{din}/all_ASCAT_segs_concat.txt",
        reference_genome="GRCh38",
        minimum_signatures=1,
        maximum_signatures=10,
        nmf_replicates=100,
        cpu=32,
    )
