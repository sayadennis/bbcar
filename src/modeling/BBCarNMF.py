import sys
import numpy as np
from sklearn.decomposition import NMF

sys.path.append('bbcar_project/src')
import BBCarModelTraining

k_list = [20, 40, 60, 80, 100, 150, 200]

# def crossval_k(X_train, y_train, X_test, y_test, k_list=k_list): # X and y should be numpy arrays. Rows of X are patients.
    
#     train_size = X_train.shape[0] # rows are patients for X
#     X = np.concatenate((X_train, X_test))
#     y = np.concatenate((y_train, y_test))
    
#     performance_dict = {}
    
#     for k in k_list:
#         A = X.T
#         nmf = NMF(n_components=k, init='nndsvd', random_state=24)
#         W = nmf.fit_transform(A)
#         F = np.dot(A.T, W)
#         F_train = F[:train_size, :]
#         F_test = F[train_size:, :]
        
#         BBCarModelTraining.record_tuning(
#             F_train, y_train, F_test, y_test, 
#             '/home/srd6051/bbcar_project/outputs/' + k + 'gene_copy_conf75_crossval_record_genethres.csv'
#         )
        
#         performance_dict[k] = [opt_mean_score, opt_params]
#     return F_test, performance_dict


def get_F(k, X_train, y_train, X_test, y_test):
    train_size = X_train.shape[0] # rows are patients for X
    X = np.concatenate((X_train, X_test))
    y = np.concatenate((y_train, y_test))

    A = X.T
    nmf = NMF(n_components=k, init='nndsvd', random_state=24)
    W = nmf.fit_transform(A)
    F = np.dot(A.T, W)
    F_train = F[:train_size, :]
    F_test = F[train_size:, :]

    return F_train, F_test, W


def coef_genes(W, genes, thres=0.01): # length of genes and rows of W should match
    if W.shape[0] != len(genes):
        return -1
    else:
        group_list = []
        for i in range(W.shape[1]): # iterate through columns of W i.e. weight vectors of each factor groups 
            coef_genes = genes[np.where(W[:,i] > thres)[0]]
            genes_final = []
            for j in range(len(coef_genes)):
                split_list = coef_genes[j].split(';') # element in list can contain multiple gene names (overlapping)
                for gene in split_list:
                    if gene in genes_final:
                        continue
                    else:
                        genes_final.append(gene)
            group_list.append([genes_final])
    return genes_final

