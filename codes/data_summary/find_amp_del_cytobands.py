import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

din = '/projects/b1122/saya/04_cleaned_cnv'
dout = '/projects/b1122/saya/05_data_summary'

# candidate_cytos = ['19q13.33', '7q22.1', '1p36.33', '1q21.1', '1q21.3'] # '4p15.33', '9p21.1' # This is from Zexian's papers
candidate_cytos = [ # these are cytobands with large LASSO coefficients 
    '12q24.33', '10q26.3', '10q11.22', '3q21.3', '2p21', # large positive coefficients
    '1q42.13', '1p32.3', '12q13.13', '15q26.3', '2p11.2'
]

cyto_copy = pd.read_csv(f'{din}/cyto_copy_conf90.csv', index_col=0)
cyto_thres_amp = pd.read_csv(f'{din}/cyto_thres_amp_conf90.csv', index_col=0)
cyto_thres_del = pd.read_csv(f'{din}/cyto_thres_del_conf90.csv', index_col=0)

case_copy = cyto_copy.iloc[['Tissue' in x for x in cyto_thres_amp.index],:]
control_copy = cyto_copy.iloc[['Control' in x for x in cyto_thres_amp.index],:]

case_amp = cyto_thres_amp.iloc[['Tissue' in x for x in cyto_thres_amp.index],:]
case_del = cyto_thres_del.iloc[['Tissue' in x for x in cyto_thres_del.index],:]
control_amp = cyto_thres_amp.iloc[['Control' in x for x in cyto_thres_amp.index],:]
control_del = cyto_thres_del.iloc[['Control' in x for x in cyto_thres_del.index],:]

for cyto in candidate_cytos:
    ## build amp/del count table
    cf = pd.DataFrame(index=['cases', 'controls'], columns=['amp only', 'del only', 'amp & del', 'none'])
    # count cases
    cf.loc['cases','amp only'] = np.sum((case_amp[cyto]>0) & (case_del[cyto]==0))
    cf.loc['cases','del only'] = np.sum((case_amp[cyto]==0) & (case_del[cyto]<0))
    cf.loc['cases', 'amp & del'] = np.sum((case_amp[cyto]>0) & (case_del[cyto]<0))
    cf.loc['cases','none'] = np.sum((case_amp[cyto]==0) & (case_del[cyto]==0))
    # count controls
    cf.loc['controls','amp only'] = np.sum((control_amp[cyto]>0) & (control_del[cyto]==0))
    cf.loc['controls','del only'] = np.sum((control_amp[cyto]==0) & (control_del[cyto]<0))
    cf.loc['controls', 'amp & del'] = np.sum((control_amp[cyto]>0) & (control_del[cyto]<0))
    cf.loc['controls','none'] = np.sum((control_amp[cyto]==0) & (control_del[cyto]==0))
    # save matrix 
    cfnorm = cf/np.transpose([cf.sum(axis=1).values])
    cfnorm.to_csv(f'{dout}/ampdelcounts_{cyto}.csv', index=True)
    ## plot histogram
    fig, ax = plt.subplots(nrows=1, ncols=1)
    ax.hist(case_copy[cyto], alpha=0.5, label='cases')
    ax.hist(control_copy[cyto], alpha=0.5, label='controls')
    ax.legend(loc='upper right')
    fig.savefig(f'{dout}/{cyto}_histogram.png')
    plt.close()

