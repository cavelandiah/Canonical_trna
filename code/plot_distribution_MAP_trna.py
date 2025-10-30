#!/usr/bin/env python3

import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

sns.set_context('talk')

experiment_map = {
    'G_i1': 'iMet_G_DMSO',
    'G_i2': 'iMet_G_NAI',
    'G_i3': 'iMet_G_DMS',
    'G_i4': 'iMet_m1G_DMSO',
    'G_i6': 'iMet_m1G_NAI',
    'G_i7': 'iMet_m1G_DMS'
}
non_vector_columns = ['ref_pos', 'ref_nt', 'align_pos',
                      'canonical_pos', 'Accessibility','Sample','Experiment_name']
samples = ['G_i1','G_i2','G_i3','G_i4','G_i6','G_i7']
all_data = []
reads_per_exp = []
thres=20
comp=sys.argv[1]
outpath_map = 'MapDistribution'
for exp in samples:
    matrix = pd.read_csv(f'{outpath_map}/concatenated_{exp}_{thres}_{comp}.csv', low_memory=False)
    matrix['Sample'] = exp
    n_reads = sum(matrix.columns.str.endswith('_vec'))
    vec_cols = [c for c in matrix.columns if c.endswith('_vec')]
    n_mutated_reads = (matrix[vec_cols] == 1).any(axis=0).sum()
    matrix['mutated_reads'] = matrix[vec_cols].sum(axis=1).astype(int)
    matrix['total_reads'] = n_reads
    matrix['total_mutated_reads'] = n_mutated_reads
    matrix['Experiment_name'] = matrix['Sample'].map(experiment_map)
    matrix[['Origin', 'Treatment']] = matrix['Experiment_name'].str.extract(r'^(iMet_[^_]+)_([^_]+)$')
    matrix['fraction_mutated'] = matrix['mutated_reads'] / matrix['total_reads']
    #matrix['fraction_mutated'] = matrix.groupby(['Origin','Treatment'])['mutated_reads'].transform(lambda x: x / x.sum())
    positions = matrix[['canonical_pos','Sample','mutated_reads','total_reads','total_mutated_reads','fraction_mutated','Origin','Treatment']]
    all_data.append(positions)

classified_df = pd.concat(all_data, ignore_index=True)
freq_df = classified_df.copy()
freq_df.to_csv(f'./MapDistribution/frequency_mutated_{comp}_{thres}.csv', index=False)
# Order treatments for consistent color mapping
treatment_order = ['DMSO', 'NAI', 'DMS']

# FacetGrid
g = sns.FacetGrid(
    freq_df,
    row='Origin',
    height=3,
    aspect=5.3,
    sharex=True,
    sharey=True,
)

# Plot mutation frequency per canonical position
g.map_dataframe(
    sns.barplot,
    x='canonical_pos',
    y='fraction_mutated',
    hue='Treatment',
    hue_order=treatment_order,
    dodge=True,
    palette='Set2'
)

g.add_legend(title='Treatment')
g.set_axis_labels("Canonical Position", "Fraction")
g.set_titles(row_template="{row_name}")
for ax in g.axes.flat:
    ax.set_ylim(0, 1)
    ax.tick_params(axis='x', rotation=90)

plt.suptitle(f"Mutation Frequency per Canonical Position ({comp}, read length={thres})", y=1.05)
plt.savefig(f"{outpath_map}/map_distribution_{comp}_{thres}_origin.pdf", bbox_inches='tight', dpi=300)
