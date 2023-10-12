import numpy as np
import pandas as pd
from inmoose.pycombat import pycombat_seq, pycombat_norm
import matplotlib.pyplot as plt 

din = '/projects/b1131/saya/bbcar/data/02b_cnv/inhouse_signatures'
dout = '/projects/b1131/saya/bbcar/data/02b_cnv/inhouse_signatures'

data = pd.read_csv(f'{din}/seglen_ampdel_category_call_counts_per_sample.csv', index_col=0)

with open('/projects/b1131/saya/bbcar/data/sample_ids_uchicago.txt', 'r') as f:
    uchicago = [int(x.strip()) for x in f.readlines()]

batches = [1 if x in uchicago else 0 for x in data.index]

data_corrected = pycombat_norm(data.values.T, batches)

pd.DataFrame(data_corrected.T, index=data.index, columns=data.columns).to_csv(f'{dout}/inhouse_cn_features_batcheffect_rm_combat.csv')

