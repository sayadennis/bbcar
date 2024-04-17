import pickle

import pandas as pd

dn = "/projects/b1122/saya/reactome_pathways"
pathway = pd.read_csv(f"{dn}/gene_to_reactome_pathway.csv")

hp = {}
for index, row in pathway.iterrows():
    pw = row["rt_pathway_id"]
    if pw not in hp:
        hp[pw] = set()
    hp[pw].add(row["hgnc_symbol"])

psz = pathway.rt_pathway_id.value_counts().to_frame()
psz.reset_index(inplace=True)
psz = psz.rename(columns={"index": "rt_pathway_id", "rt_pathway_id": "size"})
psz = psz.sort_values("size")
hpsz = {}
for index, row in psz.iterrows():
    pw = row["rt_pathway_id"]
    hpsz[pw] = row["size"]

pwids = psz["rt_pathway_id"].to_list()
contained = {}
for i in range(len(pwids) - 1):
    for j in range(len(pwids) - 1, i, -1):
        pw1, pw2 = pwids[i], pwids[j]
        if hpsz[pw1] < hpsz[pw2] and hp[pw1].issubset(hp[pw2]):
            contained[pw1] = pw2
            break

pwids = psz["rt_pathway_id"].to_list()
equaled = {}
for i in range(len(pwids) - 1):
    for j in range(len(pwids) - 1, i, -1):
        pw1, pw2 = pwids[i], pwids[j]
        if hpsz[pw1] == hpsz[pw2] and hp[pw1].issubset(hp[pw2]):
            equaled[pw1] = pw2
            break

pszf = psz.loc[
    ~psz["rt_pathway_id"].isin(list(contained.keys()) + list(equaled.keys()))
]
pszf = pszf.iloc[pszf["size"].values >= 5, :]
pszf = pszf.iloc[pszf["size"].values <= 200, :]

pszf.to_csv(f"{dn}/rt_pathway_nosubiso_ge5.csv", index=False)

# create a new filtered dictionary object and save
hp_flt = {}
for pn in pszf["rt_pathway_id"].values:
    hp_flt[pn] = hp[pn]

with open(f"{dn}/pathway_gene_dict.pkl", "wb") as f:
    pickle.dump(hp_flt, f)

## to load, run:
# with open(f'{dn}/pathway_gene_dict.pkl', 'rb') as f:
#     data = pickle.load(f)
