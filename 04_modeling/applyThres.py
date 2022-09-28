import numpy as np
import pandas as pd

def applyThres(table, thres=0.05, twotail=True):
    # table: Pandas DataFrame
    # thres: float
    bool_col = [] # list of boolean indicating the columns to keep at True
    for colname in table.columns:
        ct = np.sum(table[colname] > 0)
        rate = ct/len(table[colname])
        if twotail:
            if ((rate < thres) | (rate > (1-thres))):
                bool_col.append(False)
            else:
                bool_col.append(True)
        else:
            if (rate < thres):
                bool_col.append(False)
            else:
                bool_col.append(True)
    return table.iloc[:,bool_col]
