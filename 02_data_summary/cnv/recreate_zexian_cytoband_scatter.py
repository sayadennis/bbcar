import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt

din = '/projects/b1131/saya/bbcar/data/02b_cnv/03_gistic2_out_conf90/all'
dout = '/projects/b1131/saya/bbcar/plots/cnv'

data = pd.read_csv(f'{din}/all_lesions.conf_90.txt', sep='\t')
data['Descriptor'] = [x.strip() for x in data['Descriptor'].values]
thres = data.iloc[np.invert([x.endswith('values') for x in data['Unique Name'].values]),:data.shape[1]-1] # last column is all NaN 

genes = pd.read_csv(f'{din}/amp_genes.conf_90.txt', sep='\t')

aber_thres = {}
for aber_type in ['Amplification', 'Deletion']:
    df = thres.iloc[[x.startswith(aber_type) for x in thres['Unique Name'].values],:]
    df = pd.concat((df['Unique Name'], df.iloc[:,10:]), axis=1)
    # df = df.groupby('Descriptor').sum() # group cytobands 
    df = df.set_index('Unique Name')
    df.index.name = None
    df = (df>0).astype(int)
    aber_thres[aber_type] = df

selected_regions = {}
freq_thres = 0.75 # .75
for aber_type in ['Amplification', 'Deletion']:
    # Number of samples with amp/del in each cytoband 
    region_aber_cases = aber_thres[aber_type].iloc[:,[x.endswith('Tissue') for x in aber_thres[aber_type].columns]]
    region_aber_controls = aber_thres[aber_type].iloc[:,[x.endswith('Control') for x in aber_thres[aber_type].columns]]
    # get regions highly aberrated in cases and controls 
    highly_aber_in_cases = set(region_aber_cases.sum(axis=1).sort_values(ascending=False).iloc[region_aber_cases.sum(axis=1).sort_values(ascending=False).values>region_aber_cases.shape[1]*freq_thres].index)
    highly_aber_in_controls = set(region_aber_controls.sum(axis=1).sort_values(ascending=False).iloc[region_aber_controls.sum(axis=1).sort_values(ascending=False).values>region_aber_controls.shape[1]*freq_thres].index)
    # get regions that are aberrated only in cases 
    aber_only_in_cases = (highly_aber_in_cases - highly_aber_in_controls)
    # save to plot data
    selected_regions[aber_type] = aber_only_in_cases

plot_data = defaultdict(list)
for aber_type in ['Amplification', 'Deletion']:
    for region_name in selected_regions[aber_type]:
        corresponding_cytoband = data.iloc[data['Unique Name'].values==region_name,:].Descriptor.values[0]
        n_samples = aber_thres[aber_type].iloc[:,[x.endswith('Tissue') for x in aber_thres[aber_type].columns]].loc[region_name].sum()
        try:
            n_genes = (~pd.isnull(genes[corresponding_cytoband].iloc[3:,])).sum()
            log_q = -1* np.log10(float(genes[corresponding_cytoband].iloc[0]))
        except:
            n_genes = 0
            log_q = 0
        # 
        plot_data['amp/del'].append(aber_type)
        plot_data['region_name'].append(region_name)
        plot_data['corresponding_cytoband'].append(corresponding_cytoband)
        plot_data['n_samples'].append(n_samples)
        plot_data['n_genes'].append(n_genes)
        plot_data['log_q'].append(log_q)

#### Plot scatterplot ####
fig, ax = plt.subplots()

scatter = ax.scatter(
    plot_data['n_samples'], 
    plot_data['n_genes'], 
    s=plot_data['log_q'], 
    c=[1 if x=='Amplification' else 0 for x in plot_data['amp/del']],
    cmap='cool',
    alpha=0.5
)

for i, txt in enumerate(plot_data['corresponding_cytoband']):
    ax.annotate(txt, (plot_data['n_samples'][i], plot_data['n_genes'][i]))

# Legend 1 
# produce a legend with the unique colors from the scatter
legend1 = ax.legend(*scatter.legend_elements(),
                    loc="upper right", title="Amp/Del")
ax.add_artist(legend1)

# Legend 2 
# produce a legend with a cross section of sizes from the scatter
handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)
legend2 = ax.legend(handles, labels, loc="center right", title="-log10(q)")

ax.set_xlim(102, 118)
ax.set_xticks([105, 110, 115])
ax.set_xlabel('Number of Samples')

ax.set_ylabel('Number of Genes')

fig.savefig(f'{dout}/cnv_region_scatter_thres{freq_thres}.png')
plt.close()

###################################################
#### Plot cytobands observed in Zexian's study ####
###################################################

zexian = ['11p11.12', '10q11.21', '1q21.1', '1p36.33', '19q13.33', '9p21.1', '4p15.33', '7q22.1', '10q11.21', '7p13', '1q21.3']

thres_overlap_with_zexian = thres.iloc[[np.any([desc.startswith(x) for x in zexian]) for desc in thres['Descriptor'].values],:]

for category in ['cases', 'controls']:
    plot_data = defaultdict(list)
    for region_name in thres_overlap_with_zexian['Unique Name']:
        aber_type = region_name.split()[0]
        corresponding_cytoband = data.iloc[data['Unique Name'].values==region_name,:].Descriptor.values[0]
        if category=='cases':
            n_samples = thres.iloc[:,[x.endswith('Tissue') for x in thres.columns]].iloc[thres['Unique Name'].values==region_name].to_numpy().sum()
        elif category=='controls':
            n_samples = aber_thres[aber_type].iloc[:,[x.endswith('Control') for x in aber_thres[aber_type].columns]].loc[region_name].sum()
        try:
            n_genes = (~pd.isnull(genes[corresponding_cytoband].iloc[3:,])).sum()
            log_q = -1* np.log10(float(genes[corresponding_cytoband].iloc[0]))
        except:
            n_genes = 0
            log_q = 0
        # 
        plot_data['amp/del'].append(aber_type)
        plot_data['region_name'].append(region_name)
        plot_data['corresponding_cytoband'].append(corresponding_cytoband)
        plot_data['n_samples'].append(n_samples)
        plot_data['n_genes'].append(n_genes)
        plot_data['log_q'].append(log_q)
    # 
    #### Plot scatterplot ####
    fig, ax = plt.subplots()
    #
    scatter = ax.scatter(
        plot_data['n_samples'], 
        plot_data['n_genes'], 
        s=plot_data['log_q'], 
        c=[1 if x=='Amplification' else 0 for x in plot_data['amp/del']],
        cmap='cool',
        alpha=0.5
    )
    #
    for i, txt in enumerate(plot_data['corresponding_cytoband']):
        ax.annotate(txt, (plot_data['n_samples'][i], plot_data['n_genes'][i]), fontsize=8)
    #
    # Legend 1 
    # produce a legend with the unique colors from the scatter
    legend1 = ax.legend(*scatter.legend_elements(),
                        loc="upper right", title="Amp/Del")
    ax.add_artist(legend1)
    #
    # # Legend 2 
    # # produce a legend with a cross section of sizes from the scatter
    # handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)
    # legend2 = ax.legend(handles, labels, loc="center right", title="-log10(q)")
    #
    # ax.set_xlim(102, 118)
    # ax.set_xticks([105, 110, 115])
    ax.set_xlabel('Number of Samples')
    #
    ax.set_ylabel('Number of Genes')
    ax.set_title(f'Previously documented cytobands in {category}')
    #
    fig.savefig(f'{dout}/zexians_cytobands_in_new_dataset_scatter_{category}.png')
    pd.DataFrame(plot_data).to_csv(f'{dout}/zexians_cytobands_in_new_dataset_scatter_{category}.csv', index=False)
    plt.close()
