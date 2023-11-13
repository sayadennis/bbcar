import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import distance
from scipy.cluster import hierarchy
import seaborn as sns; sns.set_theme(color_codes=True)
from sklearn.preprocessing import StandardScaler

dout = '/projects/b1131/saya/bbcar/plots/mutation'

## Import mutational features
sbs = pd.read_csv('/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix/signature_results/SBS96/Samples.txt', sep='\t', index_col=0).T
sbs.index = [x.split('_')[0] for x in sbs.index]

dbs = pd.read_csv('/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix/signature_results/DBS78/Samples.txt', sep='\t', index_col=0).T
dbs.index = [x.split('_')[0] for x in dbs.index]

indel = pd.read_csv('/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix/signature_results/ID83/Samples.txt', sep='\t', index_col=0).T
indel.index = [x.split('_')[0] for x in indel.index]

## Import labels
label = pd.read_csv('/projects/b1131/saya/bbcar/data/clinical/bbcar_label_studyid_from_gatk_filenames.csv', index_col=0)

## Import sequencing source information
with open('/projects/b1131/saya/bbcar/data/sample_ids_uchicago.txt', 'r') as f:
    uchicago_samples = [x.strip() for x in f.readlines()]

seqsource = pd.DataFrame(np.array([1 if x in uchicago_samples else 0 for x in sbs.index]).reshape(-1,1), index=sbs.index, columns=['label'])

## Get hierarchical clustering information
data = {
    'SBS' : sbs, 'DBS' : dbs, 'INDEL' : indel
}

## Plot heatmap

for valtype in data.keys():
    mx = data[valtype]
    # scaler = StandardScaler()
    # mx = pd.DataFrame(scaler.fit_transform(mx), index=mx.index, columns=mx.columns)
    # calc row linkage 
    row_linkage = hierarchy.linkage(
        distance.pdist(np.array(mx)), method='average')
    # calc column linkage 
    col_linkage = hierarchy.linkage(
        distance.pdist(np.array(mx).T), method='average')
    # create colors for case/control labels (rows)
    rows_network_pal = sns.light_palette('lightgreen', 2)
    rows_network_lut = dict(zip([0,1], rows_network_pal))
    rows_network_colors = seqsource['label'].map(rows_network_lut)
    # create colors for signature elements (columns) 
    pos = [int(colname.startswith('+')) for colname in mx.columns]
    cols_network_pal_pos = sns.light_palette('pink', 2)
    cols_network_lut_pos = dict(zip([0,1], cols_network_pal_pos))
    cols_network_colors_pos = pd.Series(pos).map(cols_network_lut_pos)
    # plot 
    sns.clustermap(mx, row_linkage=row_linkage, col_linkage=col_linkage, 
                row_colors=rows_network_colors.values,
                col_colors=cols_network_colors_pos.values, 
                #  method="average",
                figsize=(13, 13), xticklabels=False, yticklabels=False)
    plt.savefig(f'{dout}/cluster_heatmap_{valtype}.png')
    plt.close()


