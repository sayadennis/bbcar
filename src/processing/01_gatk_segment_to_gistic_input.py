import os
import numpy as np
import pandas as pd
import itertools
import glob

din='/projects/b1122/saya/01_gatk_analyzed_segments'
dout='/projects/b1122/saya/02_gistic2_input'

combmx = pd.DataFrame(None, columns=['name', 'CONTIG', 'START', 'END', 'NUM_POINTS_COPY_RATIO', 'MEAN_LOG2_COPY_RATIO'])

for fn in glob.glob(din + '/*.csv'):
    samplen = fn.split('/')[-1].split('.')[0]
    data = pd.read_csv(fn)
    newdata = data.iloc[:,:6]
    newdata['name'] = list(itertools.repeat(samplen, len(newdata)))
    combmx = pd.concat((combmx, newdata))

combmx['MEAN_LOG2_COPY_RATIO'] = combmx['MEAN_LOG2_COPY_RATIO'] - 1

combmx.to_csv(os.path.join(dout, 'combined_gistic_input.tsv'), sep='\t', header=False, index=False)
