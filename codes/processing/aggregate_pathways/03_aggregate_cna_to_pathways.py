import numpy as np
import pandas as pd
import pickle

dn_rt = '/projects/b1122/saya/reactome_pathways'
dn_data = '/projects/b1122/saya/04_cleaned_cnv/all'
dout = '/projects/b1122/saya/08_pathway_aggregated_cnv'

#### Load the gene-to-pathway dictionary ####
with open(f'{dn_rt}/pathway_gene_dict.pkl', 'rb') as f:
    pathways = pickle.load(f)

#### Load the CNA data ####
cna = abs(pd.read_csv(f'{dn_data}/gene_thres_conf90_all.csv', index_col=0))

#### Create aggregate data ####
avg_cna = pd.DataFrame(None, index=cna.index, columns=pathways.keys())
sum_cna = pd.DataFrame(None, index=cna.index, columns=pathways.keys())
max_cna = pd.DataFrame(None, index=cna.index, columns=pathways.keys())
med_cna = pd.DataFrame(None, index=cna.index, columns=pathways.keys())

for key in pathways.keys():
    # subset cna columns by value of the dictionary 
    subcna = cna.iloc[:,[x in pathways[key] for x in cna.columns]]
    # take mean, max, median for all subjects
    avg_cna[key] = subcna.mean(axis=1)
    sum_cna[key] = subcna.sum(axis=1)
    max_cna[key] = subcna.max(axis=1)
    med_cna[key] = subcna.median(axis=1)

avg_cna.to_csv(f'{dout}/pathway_avg_cnv_studyid.csv', index=True, header=True)
sum_cna.to_csv(f'{dout}/pathway_sum_cnv_studyid.csv', index=True, header=True)
max_cna.to_csv(f'{dout}/pathway_max_cnv_studyid.csv', index=True, header=True)
med_cna.to_csv(f'{dout}/pathway_med_cnv_studyid.csv', index=True, header=True)
