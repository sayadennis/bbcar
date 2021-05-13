import os
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler

#### Load data ####
dn = '/projects/b1122/saya/additional_features'
orig_fn = 'bbcar_clinical.xlsx'

orig_data = pd.read_excel(os.path.join(dn, orig_fn), index_col=0)
sample_ids = list(pd.read_csv('bbcar/all_samples.txt', index_col=0, header=None).index)
orig_data = orig_data.loc[sample_ids]

#### Save necessary information in ML-appropriate format #### 
features = [
    'postmeno', 'benign_age', 'race_white', 'merche_age', 'been_preg', 
    'preg1age', 'rel_diag_bocancer', 'first_degree_rel_diag_bcancer'
]
data = pd.DataFrame(index=sample_ids, columns=features)

data['postmeno'] = (orig_data['Meno']=='Postmenopausal').astype(int).values
data['benign_age'] = StandardScaler().fit_transform(orig_data['Benign_Age'].values.reshape(-1,1)).reshape(-1,)
data['race_white'] = (orig_data['EthnicityCombined']=='White').astype(int).values
data['merche_age'] = np.nan_to_num(StandardScaler().fit_transform(orig_data['AgeofMerche'].values.reshape(-1,1)).reshape(-1,))
data['been_preg'] = orig_data['BeenPregnt'].values
data['preg1age'] = np.nan_to_num(StandardScaler().fit_transform(orig_data['Preg1Age'].values.reshape(-1,1)).reshape(-1,))
data['rel_diag_bocancer'] = orig_data['RelDiagBOCancer'].fillna(0).values
data['first_degree_rel_diag_bcancer'] = orig_data['FirstDegreeRelDiagBreastCancer'].fillna(0).values

data.to_csv(os.path.join(dn, 'bbcar_clinical_studyindex.csv'), header=True, index=True)

# save CSV with integer index too
data.reset_index(drop=True).to_csv(os.path.join(dn, 'bbcar_clinical_intindex.csv'), header=True, index=True)
