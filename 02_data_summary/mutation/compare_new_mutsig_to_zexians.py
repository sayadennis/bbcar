import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial import distance
from scipy.cluster import hierarchy

din = '/projects/b1131/saya/bbcar/data/02a_mutation/08_feature_matrix'
dout = '/projects/b1131/saya/bbcar/plots/mutation'

zexian = pd.read_csv('/projects/b1122/Zexian/Alignment/BBCAR_NEW/administrative/Step30Signature/Hmatrix_All.csv', index_col=0)
denovo = pd.read_csv(f'{din}/signature_results/SBS96/Suggested_Solution/SBS96_De-Novo_Solution/Signatures/SBS96_De-Novo_Signatures.txt', sep='\t', index_col=0)
cosmic = pd.read_csv(f'{din}/COSMIC_v3.3.1_SBS_GRCh38.txt', sep='\t', index_col=0)
sbs96 = pd.read_csv(f'{din}/signature_results/SBS96/Samples.txt', sep='\t', index_col=0)

def cosine_sim(a: np.array, b: np.array): 
    return np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))

#### Are any of the new COSMIC signatures similar to Zexian's O/TN? ####

fig, ax = plt.subplots(figsize=(3,4))
sim_denovo = cosine_sim(zexian.Signature_3, denovo)
ax.bar(np.arange(len(sim_denovo)), sim_denovo)
ax.axhline(y=1., linestyle='--', c='orange')
ax.set_ylim(0,1.05)
ax.set_xticks(np.arange(len(sim_denovo)))
ax.set_xticklabels(denovo.columns, rotation=60, ha='right')
ax.set_ylabel('Cosine similarity')
ax.set_title('O/TN vs. new De-Novo')
plt.tight_layout()
fig.savefig(f'{dout}/mutsig_similarity_OTN_denovo.png')
plt.close()

fig, ax = plt.subplots(figsize=(14,4))
sim_cosmic = cosine_sim(zexian.Signature_3, cosmic)
ax.bar(np.arange(len(sim_cosmic)), sim_cosmic)
ax.axhline(y=1., linestyle='--', c='orange')
ax.set_ylim(0,1.05)
ax.set_xlim(-0.5, len(sim_cosmic)-0.5)
ax.set_xticks(np.arange(len(sim_cosmic)))
ax.set_xticklabels(cosmic.columns, rotation=60, ha='right')
ax.set_ylabel('Cosine similarity')
ax.set_title('O/TN vs. COSMIC')
plt.tight_layout()
fig.savefig(f'{dout}/mutsig_similarity_OTN_cosmic.png')
plt.close()
