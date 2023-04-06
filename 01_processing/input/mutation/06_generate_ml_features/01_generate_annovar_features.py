import sys
import glob
import numpy as np
import pandas as pd
from itertools import compress

###########################################
#### Get input arguments (sample info) ####
###########################################

input_arg = sys.argv[1] # input_arg is the FULL path to VCF
source = input_arg.split('/')[-2]
sample_id = input_arg.split('/')[-1].split('.')[0].split('_')[0]
if source=='germline_only':
    pon_source = 'no'
else:
    pon_source = input_arg.split('/')[-1].split('.')[0].split('_')[1][:-3]

if sample_id.endswith('t'):
    sample_id = sample_id[:-1] # remove the "t" at the end

##########################################
#### Set input and output directories ####
##########################################

# din = '/projects/b1131/saya/bbcar/data/02a_mutation/03_annotated_variants/annovar' # included in input argument 
dout = '/projects/b1131/saya/bbcar/data/02a_mutation/04_ml_features/01_indivi_annovar_features'

#######################################
#### Set variable names to collect ####
#######################################

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

##########################
#### Collect features ####
##########################

## Define feature matrix where I can record/append data 
features = pd.DataFrame(columns=variables)

## Iterate through lines of VCF
vcf = input_arg

with open(vcf, 'r') as f:
    lines = f.readlines()

# loop through lines
for line in lines:
    if not line.startswith('#'):
        var_feature = pd.DataFrame(index=[0],columns=variables)
        # define elements of variants and their annotations 
        var_elements = line.split()
        chrom = var_elements[0]
        pos = var_elements[1]
        ref = var_elements[3]
        alt = var_elements[4]
        info = var_elements[7]
        # create variant ID 
        var_feature.loc[0,'var_id'] = f'{chrom}_{pos}_{ref}_{alt}'
        var_feature.loc[0,'source'] = source
        var_feature.loc[0,'sample_id'] = sample_id
        # get necessary information
        info = info.split(';')
        for variable in variables:
            if variable=='AF':
                af_locs = np.where([anno_item.startswith('AF=') for anno_item in info])
                afs = [x.split('=')[1] for x in np.array(info)[af_locs]]
                afs = [item for item in afs if item != '.']
                afs = [item.split(',')[0] if item=='0.5,0.5' else item for item in afs]
                if len(afs)>0:
                    af = np.mean([float(item) for item in afs])
                    var_feature.loc[0,variable] = af
            elif variable in ['var_id', 'source', 'sample_id']:
                continue
            else:
                feature_value = [item for (item, filterval) in zip(info, [anno_item.startswith(f'{variable}=') for anno_item in info]) if filterval]
                feature_value = feature_value[0].split('=')[-1]
                if feature_value.replace('.', '').replace('-', '').isnumeric():
                    var_feature.loc[0,variable] = float(feature_value)
                else:
                    var_feature.loc[0,variable] = feature_value if feature_value != '.' else None
        # record result in feature matrix 
        features = pd.concat((features,var_feature), ignore_index=True)

features.to_csv(f'{dout}/{source}_{sample_id}_{pon_source}pon_annovar_features.csv', index=False)
