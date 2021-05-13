import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

#### Load data ####
dn = '/projects/b1122/saya/additional_features'
orig_fn = 'bbcar_signature_level.xlsx'

orig_data = pd.read_excel(os.path.join(dn, orig_fn), index_col=0) # ID 467 was missing so it's imputed
sample_ids = list(pd.read_csv('bbcar/all_samples.txt', index_col=0, header=None).index)

#### Select parts of interest and standardize/scale ####
# remove study ID 467 for now since it's missing data
sample_ids.remove(467)

# only select the study IDs of interest and the four mutational signatures 
orig_data = orig_data.loc[sample_ids, ['Br_J%_1','Br_K%_3','Br_G%_30','Br_D%_MMR2']].fillna(0)

# standard scale
data = StandardScaler().fit_transform(orig_data)

# add missing subject and set all values to mean (0)
data = pd.concat(
    (pd.DataFrame(data, index=orig_data.index, columns=['Br_J%_1','Br_K%_3','Br_G%_30','Br_D%_MMR2']), 
    pd.DataFrame(np.array([[0,0,0,0]]), index=[467], columns=['Br_J%_1','Br_K%_3','Br_G%_30','Br_D%_MMR2'])), 
    axis=0
    ).sort_index()

#### Save necessary information in ML-appropriate format #### 
data.to_csv(os.path.join(dn, 'bbcar_mutsig_studyindex.csv'), header=True, index=True)

# save CSV with integer index too
data.reset_index(drop=True).to_csv(os.path.join(dn, 'bbcar_mutsig_intindex.csv'), header=True, index=True)
