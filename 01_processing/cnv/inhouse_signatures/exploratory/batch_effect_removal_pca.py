import pandas as pd
from sklearn.decomposition import PCA

din = "/projects/b1131/saya/bbcar/data/02b_cnv/inhouse_signatures"
dout = "/projects/b1131/saya/bbcar/data/02b_cnv/inhouse_signatures"

data = pd.read_csv(
    f"{din}/seglen_ampdel_category_call_counts_per_sample.csv", index_col=0
)

#############
#### PCA ####
#############

pca = PCA()
data_tf = pca.fit_transform(data.values)
data_tf[:, 0] = 0  # set PC1 to zero

data_inverse_tf = pca.inverse_transform(data_tf)

pd.DataFrame(data_inverse_tf, index=data.index, columns=data.columns).to_csv(
    f"{dout}/inhouse_sig_batcheffect_rm_pc1.csv"
)
