import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt

#### Plot scatterplot ####
fig, ax = plt.subplots(5,1, figsize=(8,18))

observed = set(['7q36.1', '19p13.2', '19q13.12', '11p15.5', '19p13.2', '19q13.42', '19q13.33', '19p13.12', '19q13.43', '1q23.1', '12p13.31', '14q32.33', '11p11.12', '17p13.1', '19p13.2', '17p12', '9p22.1'])

results = {}

for i in range(5):
    din = f'/projects/b1131/saya/bbcar/data/02b_cnv/20230126_simulate_gistic_subsample/gistic_output/subset{i}'
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
    for aber_type in ['Amplification', 'Deletion']:
        # Number of samples with amp/del in each cytoband 
        region_aber_cases = aber_thres[aber_type].iloc[:,[x.endswith('Tissue') for x in aber_thres[aber_type].columns]]
        region_aber_controls = aber_thres[aber_type].iloc[:,[x.endswith('Control') for x in aber_thres[aber_type].columns]]
        # get regions highly aberrated in cases and controls 
        freq_thres = 0
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

    scatter = ax[i].scatter(
        plot_data['n_samples'], 
        plot_data['n_genes'], 
        s=plot_data['log_q'], 
        c=[1 if x=='Amplification' else 0 for x in plot_data['amp/del']],
        cmap='cool',
        alpha=0.5
    )

    for j, txt in enumerate(plot_data['corresponding_cytoband']):
        ax[i].annotate(txt, (plot_data['n_samples'][j], plot_data['n_genes'][j]))

    # # Legend 1 
    # # produce a legend with the unique colors from the scatter
    # legend1 = ax[i].legend(*scatter.legend_elements(),
    #                     loc="upper right", title="Amp/Del")
    # ax[i].add_artist(legend1)

    # # Legend 2 
    # # produce a legend with a cross section of sizes from the scatter
    # handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)
    # legend2 = ax[i].legend(handles, labels, loc="center right", title="-log10(q)")

    # ax[i].set_xlim(102, 118)
    # ax[i].set_xticks([105, 110, 115])
    ax[i].set_xlabel('Number of Samples')
    ax[i].set_ylabel('Number of Genes')
    results[i] = pd.DataFrame(plot_data)

fig.savefig(f'{dout}/20230126_simulate_gistic_subsample.png')
plt.close()

for i in range(5):
    overlap = observed.intersection(set(results[i].corresponding_cytoband))
    print(i, len(overlap))
