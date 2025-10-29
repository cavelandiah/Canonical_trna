#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os

# Replace group labels
mapping = {
    "allowed_only": "Compatible",
    "violates_non_allowed": "Non-compatible"
}
experiment_map = {
    'G_i1': 'iMet_G_DMSO',
    'G_i2': 'iMet_G_NAI',
    'G_i3': 'iMet_G_DMS',
    'G_i4': 'iMet_m1G_DMSO',
    'G_i6': 'iMet_m1G_NAI',
    'G_i7': 'iMet_m1G_DMS'
}

threshold = sys.argv[1]
test = int(sys.argv[2])
input_folder = "../results/"
experiments = ['G_i1','G_i2','G_i3','G_i4','G_i6','G_i7']
for sample in experiments:
    if test == 1:
        path = f'{input_folder}/Test/classification_reads_{sample}_{threshold}.csv'
    else:
        path = f'{input_folder}/classification_reads_{sample}_{threshold}.csv'
    if os.path.isfile(path):
        # Smaller selection based only in vectors + reference
        summary_df = pd.read_csv(path, sep=",")
        # 1) Count groups
        group_counts = summary_df["group"].value_counts().rename_axis("group").reset_index(name="count")
        total_vecs = summary_df.shape[0]
        group_counts["proportion"] = group_counts["count"] / total_vecs
        group_counts['Sample'] = sample
        # 2) General stats (allowed vs non-allowed across *all vectors*)
        # Total counts of FP and FN
        total_FP = summary_df["FP_nonallowed_1violation"].sum()
        total_FN = summary_df["FN_allowed_0mismatch"].sum()
        total_TP = summary_df["TP_allowed_1match"].sum()
        total_TN = summary_df["TN_nonallowed_0match"].sum()

        # Total sites
        allowed_sites_total = int(summary_df["allowed_sites"].iloc[0])   # same for all vecs
        non_allowed_sites_total = int(summary_df["non_allowed_sites"].iloc[0])

        # Rates aggregated across all vectors
        overall_TPR = total_TP / (allowed_sites_total * total_vecs) if allowed_sites_total else np.nan
        overall_FPR = total_FP / (non_allowed_sites_total * total_vecs) if non_allowed_sites_total else np.nan

        general_stats = {
            "total_vecs": total_vecs,
            "n_allowed_only": (summary_df["group"] == "allowed_only").sum(),
            "n_violates_non_allowed": (summary_df["group"] == "violates_non_allowed").sum(),
            "prop_allowed_only": (summary_df["group"] == "allowed_only").mean(),
            "prop_violates_non_allowed": (summary_df["group"] == "violates_non_allowed").mean(),
            "total_TP": total_TP,
            "total_FN": total_FN,
            "total_FP": total_FP,
            "total_TN": total_TN,
            "overall_TPR_allowed": overall_TPR,
            "overall_FPR_non_allowed": overall_FPR
        }
        if test == 1:
            group_counts.to_csv(f'{input_folder}/Test/classification_{sample}_{threshold}_stats.csv', index=False)
        else:
            group_counts.to_csv(f'{input_folder}/classification_{sample}_{threshold}_stats.csv', index=False)
    else:
        continue

# PLOT
all = []
experiments = ['G_i1','G_i2','G_i3','G_i4','G_i6','G_i7']
for exp in experiments:
    if test == 1:
        fl = f'{input_folder}/Test/classification_{exp}_{threshold}_stats.csv'
    else:
        fl = f'{input_folder}/classification_{exp}_{threshold}_stats.csv'
    if os.path.isfile(fl):
        df = pd.read_csv(fl)
    else:
        continue
    all.append(df)

final_df = pd.concat(all, ignore_index=True)
final_df['Experiment_name'] = final_df['Sample'].map(experiment_map)
final_df["group"] = final_df["group"].replace(mapping)

sns.set_context('talk')
sns.barplot(final_df, x='Experiment_name', y='proportion', hue='group')
plt.xticks(rotation=45, ha='right', fontsize=10)
plt.xlabel("Experiment", fontsize=12)
plt.ylabel("Proportion", fontsize=12)
plt.title("Reference: Canonical Structure", fontsize=14)
if test == 1:
    plt.savefig(f'{input_folder}/Test/classification_compatible_mutations_{threshold}.pdf', dpi=300, bbox_inches='tight')
else:
    plt.savefig(f'{input_folder}/classification_compatible_mutations_{threshold}.pdf', dpi=300, bbox_inches='tight')

stats = final_df.groupby(['group']).agg({'proportion':'mean'})
print(stats)
