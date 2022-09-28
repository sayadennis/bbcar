import sys
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

##########################################
#### Set input and output directories ####
##########################################

din='/projects/b1042/lyglab/saya/bbcar/04_ml_features_tmp' ## CHANGE THIS ONCE b1131 IS EXPANDED 
dout='/projects/b1042/lyglab/saya/bbcar/04_ml_features_tmp' ## CHANGE THIS ONCE b1131 IS EXPANDED 
dix='/projects/b1042/lyglab/saya/bbcar/04_ml_features_tmp/somatic_pred_ix'

###################
#### Load data ####
###################

data=pd.read_csv(f'{din}/features_imputed.csv')

#####################################
#### Convert to input and target ####
#####################################

data['somatic']=(data.source=='tumor_normal').astype(int)

# nanix=data.iloc[pd.isnull(data.avsnp150).values,:].index
# data['avsnp150'].iloc[nanix]=0
# data['avsnp150'].iloc[data.avsnp150.values!=0]=1

variables=[
    'var_id',
    'AF',
    'avsnp150',
    'ExAC_ALL',
    'SIFT_score',
    'SIFT_pred',
    'Polyphen2_HDIV_score',
    'Polyphen2_HDIV_pred',
    'Polyphen2_HVAR_score',
    'Polyphen2_HVAR_pred',
    'LRT_score',
    'LRT_pred',
    'MutationTaster_score',
    'MutationTaster_pred',
    'MutationAssessor_score',
    'MutationAssessor_pred',
    'FATHMM_score',
    'FATHMM_pred',
    'MetaSVM_score',
    'MetaSVM_pred',
    'MetaLR_score',
    'MetaLR_pred',
    'VEST3_score',
    'CADD_raw',
    'CADD_phred',
    'GERP++_RS',
    'phyloP20way_mammalian',
    'phyloP100way_vertebrate',
    'SiPhy_29way_logOdds',
    'somatic'
]

Xy=data.iloc[data.source.values!='tumor_only',:][variables].drop_duplicates(ignore_index=True).set_index('var_id', drop=True)
X=Xy.iloc[:,:-1]
y=Xy.iloc[:,-1]

#######################################
#### Create train and test indices ####
#######################################

train_ix, test_ix = train_test_split(np.arange(X.shape[0]), test_size=.2, random_state=43, shuffle=True)

## Save ##
X.to_csv(f'{dout}/input.csv', index=True, header=True)
y.to_csv(f'{dout}/target.csv', index=True, header=True)

pd.DataFrame(index=train_ix).to_csv(f'{dix}/train_ix.csv', header=False)
pd.DataFrame(index=test_ix).to_csv(f'{dix}/test_ix.csv', header=False)
