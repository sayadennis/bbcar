import sys
import glob
import numpy as np
import pandas as pd

##########################################
#### Set input and output directories ####
##########################################

din='/projects/b1131/saya/bbcar/04_ml_features'
dout='/projects/b1131/saya/bbcar/04_ml_features'

############################################
#### Set empty dataframe to concatenate ####
############################################

variables=[
    'var_id',
    'source',
    'sample_id',
    'Func.refGene',
    'Gene.refGene',
    'ExonicFunc.refGene',
    'Func.knownGene',
    'Gene.knownGene',
    'GeneDetail.knownGene',
    'ExonicFunc.knownGene',
    'Func.ensGene',
    'Gene.ensGene',
    'GeneDetail.ensGene',
    'ExonicFunc.ensGene',
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
    'SiPhy_29way_logOdds'
]

features = pd.DataFrame(columns=variables)

##########################
#### Collect features ####
##########################

for fn in glob.glob(f'{din}/*_annovar_features.csv'):
    if 'bbcarpon' not in fn:
        single_sample=pd.read_csv(fn)
        features=pd.concat((features,single_sample), ignore_index=True)

# BELOW LINE MIGHT NOT BE NECESSARY if the weird '0.5,0.5' doesn't appear in AF column 
features.AF=features.AF.replace({'0.5,0.5':0.5}).astype(float)

#### Change contents of the avsnp150 column to be binary ####

nanix=features.iloc[pd.isnull(features.avsnp150).values,:].index

features['avsnp150'].iloc[nanix]=0
features['avsnp150'].iloc[features.avsnp150.values!=0]=1

## Save ##
features.to_csv(f'{dout}/annovar_features_all.csv', index=False)

##########################
#### Collect features ####
##########################

for fn in glob.glob(f'{din}/*_bbcarpon_annovar_features.csv'):
    single_sample=pd.read_csv(fn)
    features=pd.concat((features,single_sample), ignore_index=True)

# features.AF=features.AF.replace({'0.5,0.5':0.5}).astype(float)

#### Change contents of the avsnp150 column to be binary ####
nanix=features.iloc[pd.isnull(features.avsnp150).values,:].index

features['avsnp150'].iloc[nanix]=0
features['avsnp150'].iloc[features.avsnp150.values!=0]=1

## Save ##
features.to_csv(f'{dout}/annovar_features_all_bbcarpon.csv', index=False)
