# pylint: disable=wrong-import-position

import os
import sys

import numpy as np
import pandas as pd

sys.path.append("bbcar/src/modeling")
from applyThres import applyThres  # pylint: disable=import-error

din = "/projects/b1122/saya/04_cleaned_cnv"
dout = "/projects/b1122/saya/06_modified_data"

genethres = pd.read_csv(os.path.join(din, "gene_thres_conf90.csv"), index_col=0)

genethres = np.abs(genethres)
genethres.to_csv(os.path.join(dout, "gene_thres_abs.csv"))

red005 = applyThres(genethres, thres=0.05, twotail=True)
print(
    f"Shape of threshold-applied (a=0.05) matrix: {red005.shape[0]}x{red005.shape[1]}"
)

red010 = applyThres(genethres, thres=0.10, twotail=True)
print(
    f"Shape of threshold-applied (a=0.10) matrix: {red010.shape[0]}x{red010.shape[1]}"
)

red005.to_csv(os.path.join(dout, "gene_thres_thres005_abs.csv"))
red010.to_csv(os.path.join(dout, "gene_thres_thres010_abs.csv"))
