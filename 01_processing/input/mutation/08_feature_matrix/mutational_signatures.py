from SigProfilerExtractor import sigpro as sig

# Run test 
# sig.sigProfilerExtractor(input_type="vcf", output="example_output_bkp", input_data=data, 
#                          minimum_signatures=1, maximum_signatures=3, cpu=1)

# sig.sigProfilerExtractor("vcf", "example_output_results", "/projects/p30791/vcftest", reference_genome="GRCh37", 
#                          minimum_signatures=1, maximum_signatures=10, nmf_replicates=100, cpu=1)

sig.sigProfilerExtractor("vcf", "subset_results", "/projects/b1131/saya/bbcar/data/02a_mutation/07_predicted_somatic/testsubset_vcfs", 
                         reference_genome="GRCh38", minimum_signatures=1, maximum_signatures=10, nmf_replicates=100, cpu=1)

# Wait until the excecution is finished. The process may a couple of hours based on the size of the data.
# Check the current working directory for the "example_output" folder.

# sigProfilerExtractor(input_type, out_put, input_data, 
#                      reference_genome="GRCh37", opportunity_genome = "GRCh37", 
#                      context_type = "default", exome = False, 
#                      minimum_signatures=1, maximum_signatures=10, nmf_replicates=100, 
#                      resample = True, batch_size=1, cpu=-1, gpu=False, 
#                      nmf_init="random", precision= "single", matrix_normalization= "gmm", 
#                      seeds= "random", min_nmf_iterations= 10000, max_nmf_iterations=1000000, 
#                      nmf_test_conv= 10000, nmf_tolerance= 1e-15, get_all_signature_matrices= False)
