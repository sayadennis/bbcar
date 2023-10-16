import numpy as np
import pandas as pd

dout = '/projects/b1131/saya/bbcar/data'

## Load original feature matrices
mut_dn = '/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix/20230423_signature_results'
mut_fn = 'sbs_96_original_per_sample.csv'

cnv_dn = '/projects/b1131/saya/bbcar/data/02b_cnv/inhouse_signatures'
cnv_fn = 'inhouse_sig_batcheffect_rm_combat.csv'

mut = pd.read_csv(f'{mut_dn}/{mut_fn}', index_col=0)
cnv = pd.read_csv(f'{cnv_dn}/{cnv_fn}', index_col=0)

## Concatenate
combined = mut.join(cnv, how='inner')
combined.to_csv(f'{dout}/combined_mutation_cnv_orig_features.csv')

## Load signature matrices
mut_dn = '/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix'
mut_fn = 'signature_per_sample.csv'

cnv_dn = '/projects/b1131/saya/bbcar/data/02b_cnv/inhouse_signatures'
cnv_fn = 'inhouse_cnv_sig_per_sample.csv'

mut_sig = pd.read_csv(f'{mut_dn}/{mut_fn}', index_col=0)
cnv_sig = pd.read_csv(f'{cnv_dn}/{cnv_fn}', index_col=0)

## Concatenate 
combined = mut_sig.join(cnv_sig, how='inner')
combined.to_csv(f'{dout}/combined_mutation_cnv_signature_features.csv')

