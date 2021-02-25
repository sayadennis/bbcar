import sys
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt

sys.path.append("bbcar_project/src")
from LinearModelsModifier import LinearRegression

genecopy = pd.read_csv("../../data/bbcar/gene_copy_conf90.csv", index_col=0)
label = pd.read_csv("../../data/bbcar/bbcar_label.csv", index_col=0)

y = label.to_numpy().reshape(-1,1)

pvals = []
coefs = []

for colname in genecopy.columns:
    X = genecopy[colname].to_numpy().reshape(-1,1)
    lmm = LinearRegression()
    lmm.fit(X, y)
    pvals.append(lmm.p[0][0])
    coefs.append(lmm.coef_[0][0])

pval_df = pd.DataFrame(np.array([coefs, pvals]).T, index=genecopy.columns, columns=["coefs", "uncorrected p"])
pval_df_sorted = pval_df.iloc[np.argsort(pval_df["uncorrected p"].values),:]

pval_df_sorted["fdr-adjusted p"] = 0
rej, corr_p = fdrcorrection(pval_df_sorted["uncorrected p"], alpha=0.05, is_sorted=True)
pval_df_sorted["fdr-adjusted p"] = corr_p.T

pval_df_sorted.head(10).to_csv("~/bbcar_project/indv_logistic/high_pval_genes_wholegenome.csv", index=True, header=True)
