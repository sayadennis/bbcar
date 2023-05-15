import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.spatial import distance
from scipy.cluster import hierarchy

# import umap
# from sklearn.manifold import TSNE
# from sklearn.decomposition import PCA
# from sklearn.preprocessing import StandardScaler

din = '/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix'
dout = '/projects/b1131/saya/bbcar/plots/mutation'

denovo = pd.read_csv(f'{din}/signature_results/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt', sep='\t', index_col=0)
cosmic = pd.read_csv(f'{din}/COSMIC_v3.3.1_SBS_GRCh38.txt', sep='\t', index_col=0)
sbs96 = pd.read_csv(f'{din}/signature_results/SBS96/Samples.txt', sep='\t', index_col=0)
dbs78 = pd.read_csv(f'{din}/signature_results/DBS78/Samples.txt', sep='\t', index_col=0)
id83 = pd.read_csv(f'{din}/signature_results/ID83/Samples.txt', sep='\t', index_col=0)

#### Hierarchical clustering on SBS96 features ####
sbs96.columns = [x.split('_')[0] for x in sbs96]
dbs78.columns = [x.split('_')[0] for x in dbs78]
id83.columns = [x.split('_')[0] for x in id83]

col_linkage = hierarchy.linkage(distance.pdist(id83.to_numpy().T), method='average') # 'euclidean' 'seuclidean'

weird_samples = hierarchy.cut_tree(col_linkage, 2).ravel().astype(bool)
weird_samples = list(id83.loc[:,weird_samples].columns)

with open('/home/srd6051/bbcar_odd_samples_from_mutsig_10clus.txt', 'w') as f:
    for item in weird_samples:
        f.write(f'{item}\n')


