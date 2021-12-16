import os
import numpy as np
import pandas as pd
import glob

din='/projects/b1122/saya/01_gatk_analyzed_segments'
dout='/projects/b1122/saya'

sampleid_list = []
cc_list = []

for fn in glob.glob(din + '/*.csv'):
    sampleid = fn.split('/')[-1].split('.')[0].split('_')[0]
    cc = fn.split('/')[-1].split('.')[0].split('_')[1]
    sampleid_list.append(int(sampleid))
    if cc=='Tissue':
        cc_list.append(1)
    elif cc=='Control':
        cc_list.append(0)
    else:
        print('Unexpected value in case/control: {}'.format(cc))

label = pd.DataFrame(np.array(cc_list).T, index=sampleid_list, columns=['label'])
label.sort_index(axis=0, inplace=True)

label.to_csv(os.path.join(dout, 'bbcar_label_studyid.csv'), header=True, index=True)
