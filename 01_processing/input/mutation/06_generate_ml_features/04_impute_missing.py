import sys
import glob
import numpy as np
import pandas as pd

##########################################
#### Set input and output directories ####
##########################################

din='/projects/b1131/saya/bbcar/04_ml_features'
dout='/projects/b1042/lyglab/saya/bbcar/04_ml_features_tmp'

###################
#### Load data ####
###################

features=pd.read_csv(f'{din}/features_annovar_bbcarfreq.csv')

############################
#### Numerical encoding ####
############################

features.SIFT_pred=features.SIFT_pred.replace({'D':1, 'T':0})
features.Polyphen2_HDIV_pred=features.Polyphen2_HDIV_pred.replace({'B':0, 'P':1, 'D':2})
features.Polyphen2_HVAR_pred=features.Polyphen2_HVAR_pred.replace({'B':0, 'P':1, 'D':2})
features.LRT_pred=features.LRT_pred.replace({'N':0, 'U':1, 'D':2})
features.MutationTaster_pred=features.MutationTaster_pred.replace({'P':0, 'N':1, 'D':2, 'A':3})
features.MutationAssessor_pred=features.MutationAssessor_pred.replace({'N':0, 'L':1, 'M':2, 'H':3})
features.FATHMM_pred=features.FATHMM_pred.replace({'T':0, 'D':1})
features.MetaSVM_pred=features.MetaSVM_pred.replace({'T':0, 'D':1})
features.MetaLR_pred=features.MetaLR_pred.replace({'T':0, 'D':1})

####################
#### Imputation ####
####################

from sklearn.experimental import enable_iterative_imputer
from sklearn.impute import IterativeImputer

features_to_impute=[
    'AF', 'ExAC_ALL', 'SIFT_score', 'SIFT_pred', 'Polyphen2_HDIV_score',
    'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_score', 'Polyphen2_HVAR_pred',
    'LRT_score', 'LRT_pred', 'MutationTaster_score', 'MutationTaster_pred',
    'MutationAssessor_score', 'MutationAssessor_pred', 'FATHMM_score',
    'FATHMM_pred', 'MetaSVM_score', 'MetaSVM_pred', 'MetaLR_score',
    'MetaLR_pred', 'VEST3_score', 'CADD_raw', 'CADD_phred', 'GERP++_RS',
    'phyloP20way_mammalian', 'phyloP100way_vertebrate',
    'SiPhy_29way_logOdds'
]

imp_median=IterativeImputer(random_state=0, initial_strategy='median')
imp_features=imp_median.fit_transform(features[features_to_impute])

features[features_to_impute]=imp_features

## Save ##
features.to_csv(f'{dout}/features_imputed.csv', index=False, header=True)
