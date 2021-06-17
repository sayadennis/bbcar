import sys
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
from statsmodels.stats.multitest import fdrcorrection
import matplotlib.pyplot as plt

genethres = np.absolute(pd.read_csv("../../data/bbcar/gene_thres_conf90.csv", index_col=0))
regthres = np.absolute(pd.read_csv("../../data/bbcar/reg_thres_conf90.csv", index_col=0).T)
label = pd.read_csv("../../data/bbcar/bbcar_label.csv", index_col=0)

y = label.to_numpy().reshape(-1,1)


# First I will test the gene-level threshold features 

gene_pvals = []
gene_chis = []

genethres_pos = genethres.iloc[y==1,:]
genethres_neg = genethres.iloc[y==0,:]

for colname in genethres.columns:
    mx = np.array(
        [[np.sum([x != 0 for x in genethres_pos[colname]]), np.sum([x != 0 for x in genethres_neg[colname]])],
        [np.sum([x == 0 for x in genethres_pos[colname]]), np.sum([x == 0 for x in genethres_neg[colname]])]]
    )
    chival, pval, _, _ = chi2_contingency(mx)
    gene_pvals.append(pval)
    gene_chis.append(chival)

pval_df = pd.DataFrame(np.array([gene_chis, gene_pvals]).T, index=genethres.columns, columns=["chi value", "uncorrected p"])
pval_df_sorted = pval_df.iloc[np.argsort(pval_df["uncorrected p"].values),:]

pval_df_sorted["fdr-adjusted p"] = 0
rej, corr_p = fdrcorrection(pval_df_sorted["uncorrected p"], alpha=0.05, is_sorted=True)
pval_df_sorted["fdr-adjusted p"] = corr_p.T

pval_df_sorted.head(10).to_csv("~/bbcar_project/thres_chitest/chi2_high_pval_genethres.csv", index=True, header=True)

# Next I will test the region-level threshold features

reg_pvals = []
reg_chis = []

regthres_pos = regthres.iloc[y==1,:]
regthres_neg = regthres.iloc[y==0,:]

for colname in regthres.columns:
    mx = np.array(
        [[np.sum([x != 0 for x in regthres_pos[colname]]), np.sum([x != 0 for x in regthres_neg[colname]])],
        [np.sum([x == 0 for x in regthres_pos[colname]]), np.sum([x == 0 for x in regthres_neg[colname]])]]
    )
    chival, pval, _, _ = chi2_contingency(mx)
    reg_pvals.append(pval)
    reg_chis.append(chival)

pval_df = pd.DataFrame(np.array([reg_chis, reg_pvals]).T, index=regthres.columns, columns=["chi value", "uncorrected p"])
pval_df_sorted = pval_df.iloc[np.argsort(pval_df["uncorrected p"].values),:]

pval_df_sorted["fdr-adjusted p"] = 0
rej, corr_p = fdrcorrection(pval_df_sorted["uncorrected p"], alpha=0.05, is_sorted=True)
pval_df_sorted["fdr-adjusted p"] = corr_p.T

pval_df_sorted.head(10).to_csv("~/bbcar_project/thres_chitest/chi2_high_pval_regthres.csv", index=True, header=True)
