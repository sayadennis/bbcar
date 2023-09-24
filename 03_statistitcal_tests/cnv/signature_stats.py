import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

################################
#### Get signature exposure ####
################################

fin = '/projects/b1131/saya/bbcar/data/02b_cnv/signatures/04_signatures/CNV48/Samples.txt'

cn_features = pd.read_csv(fin, sep='\t')

cosmic_features = pd.read_csv('/projects/b1131/saya/bbcar/data/02b_cnv/signatures/04_signatures/CNV48/Suggested_Solution/COSMIC_CNV48_Decomposed_Solution/Signatures/COSMIC_CNV48_Signatures.txt', sep='\t')

exposures = pd.DataFrame(np.dot(cn_features.iloc[:,1:].T, cosmic_features.iloc[:,1:]), index=cn_features.columns[1:], columns=cosmic_features.columns[1:])

###############################
#### Get case/control info ####
###############################

outcome = pd.read_csv('/projects/b1131/saya/bbcar/data/clinical/bbcar_label_studyid_from_gatk_filenames.csv', index_col=0)
case_ids = [str(ix) for ix in outcome.index if outcome.loc[ix,'label']==1]

#################################
#### Compare exposure levels ####
#################################

for sig in exposures.columns:
    case_exposures = exposures.iloc[[x in case_ids for x in exposures.index],:][sig].values.ravel()
    control_exposures = exposures.iloc[[x not in case_ids for x in exposures.index],:][sig].values.ravel()
    t, p = ttest_ind(case_exposures, control_exposures)
    print(f'Signature {sig}: t={t:.2f}, p={p:.4f}')
