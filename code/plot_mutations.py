#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt


all_data = []
thres = 75
experiments = ['G_i1','G_i2','G_i3','G_i4','G_i6','G_i7']
for f in experiments:
    df = pd.read_csv(f'Results/distance_{f}_{thres}.csv')
    all_data.append(df)
cdf = pd.concat(all_data, ignore_index=True)

experiment_map = {
    'G_i1': 'iMet_G_DMSO',
    'G_i2': 'iMet_G_NAI',
    'G_i3': 'iMet_G_DMS',
    'G_i4': 'iMet_m1G_DMSO',
    'G_i6': 'iMet_m1G_NAI',
    'G_i7': 'iMet_m1G_DMS'
}
cdf['Experiment_name'] = cdf['Experiment'].map(experiment_map)
cdf['Origin'] = cdf['Experiment_name'].str.split('_').str[:2].str.join('_')
cdf['Treatment'] = cdf['Experiment_name'].str.split('_').str[-1:].str.join('_')
cdf['mutations'] = cdf['mutations'].astype(int)

cdf_original = cdf.copy()
cdf = cdf.dropna()

sns.set(style="whitegrid", context="talk")

cdf_original.groupby(['Experiment_name']).agg({'mutations':'mean',
                                     'mutations': 'min'})
# Mutation dist
f = sns.displot(
    data=cdf_original,
    x='mutations',
    col ='Origin',
    hue='Treatment',
    kind="hist",
    multiple="layer",
    element='step',
    stat='density',
    discrete=True # Add this since should be a category in the histplot
)
# Add log scale
#for ax in f.axes.flat:
#    ax.set_yscale("log")

f._legend.set_bbox_to_anchor((1, 0.5))  # Move to the right
f._legend.set_title("Treatment")  # Optional: rename title

#plt.xticks(rotation=0, ha='right', fontsize=10)
#plt.xlabel("Mutation Number", fontsize=12)
#plt.ylabel("Similarity to canonical str.", fontsize=12)
plt.savefig(f'mutations_per_read_{thres}.pdf', dpi=300, bbox_inches='tight')
