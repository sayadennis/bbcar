import os
import numpy as np
import pandas as pd

# Load data 
label = pd.read_csv('/projects/b1122/saya/bbcar_label_studyid.csv', index_col=0)
pon = pd.read_csv('/home/srd6051/bbcar/sample_ids_pon.txt', index_col=0)

# Get labels of only PON 
pon_label = label.iloc[[i in pon.index for i in label.index],:]
# Get labels of only non-PON
nonpon_label = label.iloc[[i not in pon.index for i in label.index],:]

# Record the case/control and PON/non-PON sample composition 
composition = pd.DataFrame(0, index=['case', 'control'], columns=['tissue', 'tissue & germline'])

composition.loc['case', 'tissue & germline'] = np.sum(pon_label.values==1)
composition.loc['control', 'tissue & germline'] = np.sum(pon_label.values==0)
composition.loc['case', 'tissue'] = np.sum(nonpon_label.values==1)
composition.loc['control', 'tissue'] = np.sum(nonpon_label.values==0)

print(composition)

# print('Case : Control = %s : %s' % (np.sum(label.values), np.sum(label.values==0)))
