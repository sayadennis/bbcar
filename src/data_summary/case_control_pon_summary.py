import os
import numpy as np
import pandas as pd

label = pd.read_csv('/projects/b1122/saya/bbcar_label_studyidindex.csv', index_col=0)
pon = pd.read_csv('/home/srd6051/bbcar/sample_ids_pon.txt', index_col=0)

print('Case : Control = %s : %s' % (np.sum(label.values), np.sum(label.values==0)))
