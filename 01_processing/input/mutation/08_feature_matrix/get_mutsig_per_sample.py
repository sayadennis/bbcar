import numpy as np
import pandas as pd

din = '/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix/signature_results/SBS96'
dout = ''

samples_sbs96 = pd.read_csv(f'{din}/Samples.txt', sep='\t', index_col=0)
s5_solution = pd.read_csv(f'{din}/All_Solutions/SBS96_5_Signatures/Signatures/SBS96_S5_Signatures.txt', sep='\t', index_col=0)

sigs_per_sample = np.dot(samples_sbs96, s5_solution) # realized that I had duplicate VCFs for matched samples - worknig on this now

