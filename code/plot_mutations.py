#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# --- Count mutations per Experiment_name ---
def classify_mut(x):
    if x == 0:
        return "0"
    elif x == 1:
        return "1"
    else:
        return ">1"

# Input
thres = sys.argv[1]
test= sys.argv[2]

if test == 1:
    input_folder='../results/Test'
    # Output
    output_folder="../results/Test"
else:
    input_folder='../results'
    # Output
    output_folder="../results"

all_data = []
experiments = ['G_i1','G_i2','G_i3','G_i4','G_i6','G_i7']
for f in experiments:
    file_stats = f'{input_folder}/distance_{f}_{thres}.csv'
    if os.path.isfile(file_stats):
        df = pd.read_csv(file_stats)
        all_data.append(df)
    else:
        continue
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
plt.savefig(f'{output_folder}/mutations_per_read_{thres}.pdf', dpi=300, bbox_inches='tight')

## Table:

# Add a categorical column
cdf_original["mut_class"] = cdf_original["mutations"].apply(classify_mut)

# Count by Experiment_name and mutation class
counts = (
    cdf_original
    .groupby(["Experiment_name", "mut_class"])
    .size()
    .unstack(fill_value=0)
    .reset_index()
)

# Add total reads per experiment
counts["total_reads"] = counts[["0", "1", ">1"]].sum(axis=1)

# Optional: reorder columns
counts = counts[["Experiment_name", "0", "1", ">1", "total_reads"]]

print("\n=== Read counts per experiment ===")
print(counts.to_string(index=False))

# Optionally, write table to CSV
counts.to_csv(f"{output_folder}/mutation_summary_{thres}.csv", index=False)

