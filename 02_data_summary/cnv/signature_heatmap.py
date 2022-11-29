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

din = '/projects/b1131/saya/bbcar/data/02b_cnv/signatures'
dout = '/projects/b1131/saya/bbcar/plots/cnv'

data = {
    'counts' : {
        'length' : pd.read_csv(f'{din}/seglen_category_call_counts_per_sample.csv', index_col=0, header=[0,1]),
        'levels' : pd.read_csv(f'{din}/seglen_ampdel_category_call_counts_per_sample.csv', index_col=0, header=[0,1,2]),
    },
    'ratios' : {
        'length' : pd.read_csv(f'{din}/seglen_category_call_ratios_per_sample.csv', index_col=0, header=[0,1]),
        'levels' : pd.read_csv(f'{din}/seglen_ampdel_category_call_ratios_per_sample.csv', index_col=0, header=[0,1,2]),
    },
}

for valtype in ['counts', 'ratios']:
    for cattype in ['length', 'levels']:
        mx = data[valtype][cattype]
        scaler = StandardScaler()
        mx_scaled = scaler.fit_transform(mx)
        # calc row linkage 
        row_linkage = hierarchy.linkage(
            distance.pdist(np.array(mx_scaled)), method='average')
        # calc column linkage 
        col_linkage = hierarchy.linkage(
            distance.pdist(np.array(mx_scaled).T), method='average')
        # create colors for case/control labels (rows)
        labels = np.array(['Tissue' in samplename for samplename in mx.index], dtype=float)
        rows_network_pal = sns.light_palette('orange', 2)
        rows_network_lut = dict(zip([0,1], rows_network_pal))
        rows_network_colors = pd.Series(labels).map(rows_network_lut)
        # create colors for signature elements (columns) 
        #### FIX HERE - NEED MULTIINDEX NAMES TO CORRECTLY GET COLORS ####
        cols_network_pal1 = sns.light_palette('red', len(mx.iloc[:,mx.columns.get_level_values(1)=='+'].columns.get_level_values(1).unique()))
        cols_network_pal2 = sns.light_palette('blue', len(mx.iloc[:,mx.columns.get_level_values(1)=='-'].columns.get_level_values(1).unique()))
        cols_network_lut1 = dict(zip(mx.iloc[:,mx.columns.get_level_values(1)=='+'].columns.get_level_values(0).unique(), cols_network_pal1))
        cols_network_lut2 = dict(zip(mx.iloc[:,mx.columns.get_level_values(1)=='-'].columns.get_level_values(0).unique(), cols_network_pal2))
        cols_network_colors = pd.Series(labels, dtype=object)
        ##################################################################
        for i in cols_network_colors.index:
            if mx.iloc[:,i].name[1]=='+':
                cols_network_colors.loc[i] = cols_network_lut1[mx.iloc[:,i].name[0]]
            elif mx.iloc[:,i].name[1]=='-':
                cols_network_colors.loc[i] = cols_network_lut2[mx.iloc[:,i].name[0]]
        # plot 
        sns.clustermap(mx_scaled, row_linkage=row_linkage, col_linkage=col_linkage, 
                    row_colors=rows_network_colors.values,
                    col_colors=cols_network_colors.values, 
                    #  method="average",
                    figsize=(13, 13), xticklabels=False, yticklabels=False)
        plt.savefig(f'{dout}/cluster_heatmap_{valtype}_{cattype}.png')
        plt.close()

