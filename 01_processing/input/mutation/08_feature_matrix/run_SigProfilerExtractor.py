from SigProfilerExtractor import sigpro as sig

if __name__ == "__main__":
    sig.sigProfilerExtractor(input_type="vcf", 
                             output="/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix/signature_results", 
                             input_data="/projects/b1131/saya/bbcar/data/02a_mutation/07_predicted_somatic/vcfs", 
                             reference_genome="GRCh38", minimum_signatures=1, 
                             maximum_signatures=10, nmf_replicates=100, cpu=32)
