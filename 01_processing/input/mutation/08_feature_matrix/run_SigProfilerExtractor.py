from SigProfilerExtractor import sigpro as sig

din = '/projects/b1131/saya/bbcar/data/02a_mutation/07_predicted_somatic/vcfs/bbcar'
dout = '/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix/20230423_signature_results'

if __name__ == "__main__":
    sig.sigProfilerExtractor(input_type='vcf', 
                             output=dout,
                             input_data=din,
                            #  input_data="/projects/b1131/saya/bbcar/data/02a_mutation/07_predicted_somatic/vcfs_matched", 
                             reference_genome='GRCh38', minimum_signatures=1, 
                             maximum_signatures=10, nmf_replicates=100, cpu=32)
