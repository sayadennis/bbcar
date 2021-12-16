import os
import sys
import numpy as np
import pandas as pd
import glob

sys.path.append('bbcar/src/processing')
from process_functions import lesions_to_structured_regions
from process_functions import generate_ovelap_mx
from process_functions import transform_gatk_to_reg_matrix

d01='/projects/b1122/saya/01_gatk_analyzed_segments'
d03='/projects/b1122/saya/03_gistic2_out_conf90'
d04='/projects/b1122/saya/04_cleaned_cnv'
dix='/projects/b1122/saya/indices'

## Get test-set index
all_ix = pd.read_csv(f'{dix}/index_dictionary.csv', index_col=0)
all_studyid = all_ix['studyid']

test_ix = pd.read_csv(f'{dix}/test_ix.csv', header=None)[0] # this is index 0-200 that I use in the "reindex" step 
test_studyid = pd.read_csv(f'{dix}/test_studyid.csv', header=None)[0] # this is the study ID that matches the GATK file names

## Get cytoband coordinates and re-format so that we can use it in "transform_gatk_to_reg_matrix" function later
cyto_coords = pd.read_csv('/projects/b1122/saya/ref_cytobands/cytoband_hg38.tsv', sep='\t')
cyto_coords = cyto_coords.iloc[~pd.isnull(cyto_coords['name'].values),:].reset_index(drop=True) # only select cytobands with names (characterized)
cyto_regions = pd.DataFrame({
    'chrom' : cyto_coords['#chrom'],
    'start' : cyto_coords['chromStart'],
    'end' : cyto_coords['chromEnd']
})
cyto_regions.index = [cyto_coords['#chrom'].iloc[i].split('chr')[1] + cyto_coords['name'].iloc[i] for i in cyto_coords.index]

## Get gene coordinates and re-format so that we can use it in "transform_gatk_to_reg_matrix" function later
gene_coords = pd.read_csv('/projects/b1122/saya/ref_gene_annotations/GRCh38_gencode.v27.refFlat.txt', sep='\t', header=None)
chroms = list(cyto_regions['chrom'].unique())
gene_coords = gene_coords.iloc[[x in chroms for x in gene_coords[2].values],:].reset_index(drop=True) # only select normal chromosomes
gene_coords = gene_coords[~gene_coords[0].duplicated(keep='first')].reset_index(drop=True) # only keep one instance of each gene
gene_regions = pd.DataFrame({
    'chrom' : gene_coords[2],
    'start' : gene_coords[4],
    'end' : gene_coords[5]
})
gene_regions.index = [x for x in gene_coords[0].values]


###############################################
#### All cases and all controls separately ####
###############################################

case_regions = lesions_to_structured_regions(f'{d03}/all_cases/all_lesions.conf_90.txt')
control_regions = lesions_to_structured_regions(f'{d03}/all_controls/all_lesions.conf_90.txt')

ol_mx = generate_ovelap_mx(case_regions, control_regions)

ol_thres = 1e+2 # minimum amount of overlap between regions to consider "overlapping" 
case_unique_regions = list(ol_mx.index[ol_mx.sum(axis=1).values<ol_thres])
control_unique_regions = list(ol_mx.columns[ol_mx.sum(axis=0).values<ol_thres])

# create a reg_table with these unique regions
reg_table = pd.concat(
    (case_regions.iloc[[x in case_unique_regions for x in case_regions.index]], 
    control_regions.iloc[[x in control_unique_regions for x in control_regions.index]]),
    axis=0
)
# reindex with case/control info to make indices unique
reg_table.index = ['Case ' + x for x in case_unique_regions] + ['Control ' + x for x in control_unique_regions]

allsample_mx = transform_gatk_to_reg_matrix(reg_table=reg_table, gatk_dir=d01)
allsample_mx.to_csv(f'{d04}/all_ccsep/reg_copy_unique.csv', header=True, index=True)


#########################################################
#### Training cases and training controls separately ####
#########################################################

## Cytoband-level features - record the amp/del with largest absolute value within each cytoband 
train_cytocopy = pd.read_csv(f'{d04}/train_ccsep/cyto_copy_conf90_ccsep.csv', index_col=0) # read training matrix 
aber_cyto_regions = cyto_regions.iloc[[x in list(train_cytocopy.columns) for x in cyto_regions.index],:]

test_cytocopy = transform_gatk_to_reg_matrix(aber_cyto_regions, d01, sampleset=list(test_studyid))

ccsep_cytocopy = pd.concat((train_cytocopy, test_cytocopy), axis=0).sort_index()
ccsep_cytocopy.to_csv(f'{d04}/train_ccsep/cyto_copy.csv', header=True, index=True)

## Gene-level features - if a sample has amp/del segment overlapping a gene, enter the log copy number 
# Get gene feature list from cleaned training set
train_genecopy = pd.read_csv(f'{d04}/train_ccsep/gene_copy_conf90_ccsep.csv', index_col=0) # read training matrix 
gene_regions = gene_regions.iloc[[x in list(train_genecopy.columns) for x in gene_regions.index],:]

test_genecopy = transform_gatk_to_reg_matrix(gene_coords, d01, sampleset=list(test_studyid))

ccsep_genecopy = pd.concat((train_genecopy, test_genecopy), axis=0).sort_index()
ccsep_genecopy.to_csv(f'{d04}/train_ccsep/gene_copy.csv', header=True, index=True)

## Region-level features - get case-unique and control-unique regions, record CN value of overlapping segments
# read case and control regions
case_regions = lesions_to_structured_regions(f'{d03}/train_cases/all_lesions.conf_90.txt')
control_regions = lesions_to_structured_regions(f'{d03}/train_controls/all_lesions.conf_90.txt')
# compute overlaps
ol_mx = generate_ovelap_mx(case_regions, control_regions)
# find number of overlapping regions for case & control respectively (unnecessary)
olsize_control = np.sum(ol_mx.sum(axis=0).values>ol_thres)
olsize_case = np.sum(ol_mx.sum(axis=1).values>ol_thres)
# get a list of case-unique and control-unique region names 
case_unique_regions = list(ol_mx.index[ol_mx.sum(axis=1).values<ol_thres])
control_unique_regions = list(ol_mx.columns[ol_mx.sum(axis=0).values<ol_thres])
# create a reg_table with these unique regions
reg_table = pd.concat(
    (case_regions.iloc[[x in case_unique_regions for x in case_regions.index]], 
    control_regions.iloc[[x in control_unique_regions for x in control_regions.index]]),
    axis=0
)
# reindex with case/control info to make row-indices unique
reg_table.index = ['Case ' + x for x in case_unique_regions] + ['Control ' + x for x in control_unique_regions]
# save
allsample_mx = transform_gatk_to_reg_matrix(reg_table=reg_table, gatk_dir=d01)
allsample_mx.to_csv(f'{d04}/train_ccsep/reg_cc_unique_copy_gatk_all.csv', header=True, index=True)


###################################################################################
#### (Training cases + test all) and (training controls + test all) separately ####
###################################################################################

## Cytoband-level features 
cytocopy = pd.read_csv(f'{d04}/train_ccsep_test_all/cyto_copy_conf90_ccsep.csv', index_col=0)
cytocopy_transformed = cytocopy.iloc[[int(x.split('_')[0]) not in test_studyid for x in cytocopy.index], :]
cytocopy_test = cytocopy.iloc[[int(x.split('_')[0]) in test_studyid for x in cytocopy.index], :]
for studyix in cytocopy_test.index.unique():
    subset = cytocopy_test.iloc[cytocopy_test.index==studyix,:]
    cytocopy_transformed.loc[studyix,:] = list(subset.mean())

cytocopy_transformed.to_csv(f'{d04}/train_ccsep_test_all/cyto_copy.csv', header=True, index=True)

## Gene-level features
genecopy = pd.read_csv(f'{d04}/train_ccsep_test_all/gene_copy_conf90_ccsep.csv', index_col=0)
genecopy_transformed = genecopy.iloc[[int(x.split('_')[0]) not in test_studyid for x in genecopy.index], :]
genecopy_test = genecopy.iloc[[int(x.split('_')[0]) in test_studyid for x in genecopy.index], :]
for studyix in genecopy_test.index.unique():
    subset = genecopy_test.iloc[genecopy_test.index==studyix,:]
    genecopy_transformed.loc[studyix,:] = list(subset.mean())

genecopy_transformed.to_csv(f'{d04}/train_ccsep_test_all/gene_copy.csv', header=True, index=True)

## Region-level features
# read case and control regions
case_regions = lesions_to_structured_regions(f'{d03}/train_cases_test_all/all_lesions.conf_90.txt')
control_regions = lesions_to_structured_regions(f'{d03}/train_controls_test_all/all_lesions.conf_90.txt')
# compute overlaps
ol_mx = generate_ovelap_mx(case_regions, control_regions)
# find number of overlapping regions for case & control respectively (unnecessary)
olsize_control = np.sum(ol_mx.sum(axis=0).values!=0)
olsize_case = np.sum(ol_mx.sum(axis=1).values!=0)
# get a list of case-unique and control-unique region names 
case_unique_regions = list(ol_mx.index[ol_mx.sum(axis=1).values==0])
control_unique_regions = list(ol_mx.columns[ol_mx.sum(axis=0).values==0])
# create a reg_table with these unique regions
reg_table = pd.concat(
    (case_regions.iloc[[x in case_unique_regions for x in case_regions.index]], 
    control_regions.iloc[[x in control_unique_regions for x in control_regions.index]]),
    axis=0
)
# reindex with case/control info to make indices unique
reg_table.index = ['Case ' + x for x in case_unique_regions] + ['Control ' + x for x in control_unique_regions]
# save
allsample_mx = transform_gatk_to_reg_matrix(reg_table=reg_table, gatk_dir=d01)
allsample_mx.to_csv(f'{d04}/train_ccsep_test_all/reg_copy_unique.csv', header=True, index=True)
