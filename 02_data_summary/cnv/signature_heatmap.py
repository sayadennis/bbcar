import numpy as np
import pandas as pd
from scipy.spatial import distance
from scipy.cluster import hierarchy
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_theme(color_codes=True)
from sklearn.preprocessing import StandardScaler

###################
#### Load data ####
###################

din = '/projects/b1131/saya/bbcar/data/02b_cnv'
dout = '/projects/b1131/saya/bbcar/plots/cnv'

data = {
    'counts' : pd.read_csv(f'{din}/signatures/seglen_ampdel_category_call_counts_per_sample.csv', index_col=0, header=0),
    'ratios' : pd.read_csv(f'{din}/signatures/seglen_ampdel_category_call_ratios_per_sample.csv', index_col=0, header=0),
    'counts_pca' : pd.read_csv(f'{din}/inhouse_signatures/inhouse_sig_batcheffect_rm_pc1.csv', index_col=0, header=0),
    'counts_combat' : pd.read_csv(f'{din}/inhouse_signatures/inhouse_sig_batcheffect_rm_combat.csv', index_col=0, header=0),
}

# labels = pd.read_csv(f'/projects/b1131/saya/bbcar/data/clinical/bbcar_label_studyid_from_gatk_filenames.csv', index_col=0)
# labels = labels.loc[data['counts'].index]

## labels to use: U Chicago vs. Indiana sequencing (observing batch effect)
with open('/projects/b1131/saya/bbcar/data/sample_ids_uchicago.txt', 'r') as f:
    uchicago_samples = [int(x.strip()) for x in f.readlines()]

labels = pd.DataFrame(np.array([1 if x in uchicago_samples else 0 for x in data['counts'].index]).reshape(-1,1), index=data['counts'].index, columns=['label'])

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
    rows_network_colors = labels['label'].map(rows_network_lut)
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
    plt.savefig(f'{dout}/cluster_heatmap_{valtype}_notstandardscaled.png')
    plt.close()

valtype = 'counts'

mx = data[valtype]
row_linkage = hierarchy.linkage(
    distance.pdist(np.array(mx)), method='average')

weird_samples = hierarchy.cut_tree(row_linkage, 2).ravel().astype(bool)
weird_samples = list(mx.iloc[weird_samples,:].index)

with open('/home/srd6051/bbcar_odd_samples.txt', 'w') as f:
    for item in weird_samples:
        f.write(f'{item}\n')
