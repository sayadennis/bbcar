# import os
# import sys
# import numpy as np
# import pandas as pd
from SigProfilerExtractor import sigpro as sig

# to get input from vcf files
path_to_example_folder_containing_vcf_files = sig.importdata("vcf")

# you can put the path to your folder containing the vcf samples
data = path_to_example_folder_containing_vcf_files

# Run test 
sig.sigProfilerExtractor(input_type="vcf", output="example_output", input_data=data, 
                         minimum_signatures=1, maximum_signatures=3)

# Wait until the excecution is finished. The process may a couple of hours based on the size of the data.
# Check the current working directory for the "example_output" folder.

# sigProfilerExtractor(input_type, output, input_data, 
#                      reference_genome="GRCh37", opportunity_genome = "GRCh37", 
#                      context_type = "default", exome = False, 
#                      minimum_signatures=1, maximum_signatures=10, nmf_replicates=100, 
#                      resample = True, batch_size=1, cpu=-1, gpu=False, nmf_init="random", 
#                      precision= "single", matrix_normalization= "gmm", seeds= "random", 
#                      min_nmf_iterations= 10000, max_nmf_iterations=1000000, 
#                      nmf_test_conv= 10000, nmf_tolerance= 1e-15, 
#                      get_all_signature_matrices= False)
