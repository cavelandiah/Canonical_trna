#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os

sns.set_context('talk')
#sns.set_style('whitegrid')

# ---------------------------
# Input arguments and setup
# ---------------------------
threshold = sys.argv[1]
test = int(sys.argv[2])

input_folder = "../results/"
subfolder = "Test" if test == 1 else ""
prefix = f"{input_folder}/{subfolder}" if subfolder else input_folder

experiments = ['G_i1','G_i2','G_i3','G_i4','G_i6','G_i7']

mapping = {
    "Compatible": "Compatible",
    "Noncompatible": "Non-compatible",
    "Non-informative":'Non-informative'
}
experiment_map = {
    'G_i1': 'iMet_G_DMSO',
    'G_i2': 'iMet_G_NAI',
    'G_i3': 'iMet_G_DMS',
    'G_i4': 'iMet_m1G_DMSO',
    'G_i6': 'iMet_m1G_NAI',
    'G_i7': 'iMet_m1G_DMS'
}

# ---------------------------
# Collect summaries
# ---------------------------
all_stats = []
for sample in experiments:
    path = f'{prefix}/classification_reads_{sample}_{threshold}.csv'
    if not os.path.isfile(path):
        continue

    df = pd.read_csv(path)
    total_vecs = len(df)

    # --- Count groups ---
    group_counts = (
        df["group"].value_counts()
        .rename_axis("group")
        .reset_index(name="count")
    )
    group_counts["proportion"] = group_counts["count"] / total_vecs
    group_counts['Sample'] = sample

    # --- General TP/FP summary ---
    total_FP = df["FP"].sum()
    total_FN = df["FN"].sum()
    total_TP = df["TP"].sum()
    total_TN = df["TN"].sum()

    allowed_sites_total = int(df["allowed_sites"].iloc[0])
    non_allowed_sites_total = int(df["non_allowed_sites"].iloc[0])

    overall_TPR = total_TP / (allowed_sites_total * total_vecs) if allowed_sites_total else np.nan
    overall_FPR = total_FP / (non_allowed_sites_total * total_vecs) if non_allowed_sites_total else np.nan

    # --- Append to master list ---
    group_counts["TPR_allowed"] = overall_TPR
    group_counts["FPR_non_allowed"] = overall_FPR
    all_stats.append(group_counts)

# ---------------------------
# Combine results
# ---------------------------
if not all_stats:
    sys.exit("No valid result files found.")

final_df = pd.concat(all_stats, ignore_index=True)
final_df["group"] = final_df["group"].replace(mapping)
final_df["Experiment_name"] = final_df["Sample"].map(experiment_map)

# Save combined CSV
out_csv = f'{prefix}/classification_summary_{threshold}.csv'
final_df.to_csv(out_csv, index=False)

# ---------------------------
# Summary table (mean ± SD)
# ---------------------------
summary_stats = (
    final_df.groupby(["group"])["proportion"]
    .agg(['mean', 'std'])
    .rename(columns={'mean': 'mean_proportion', 'std': 'sd_proportion'})
)
print("\n=== Summary of group proportions ===")
print(summary_stats.round(3))
print("\nCombined CSV saved to:", out_csv)

# ---------------------------
# Plot
# ---------------------------
plt.figure(figsize=(8, 5))
ax = sns.barplot(
    data=final_df,
    x='Experiment_name',
    y='proportion',
    hue='group',
    hue_order=['Non-informative','Non-compatible','Compatible']
    #errorbar='sd'
)
plt.xticks(rotation=45, ha='right', fontsize=10)
plt.xlabel("Experiment", fontsize=12)
plt.ylabel("Proportion", fontsize=12)
plt.title("Reference: Canonical Structure", fontsize=14)
plt.legend(title='Group')
ax.legend(
    title='Group',
    bbox_to_anchor=(1.05, 1),  # x=1.05 pushes it to the right
    loc='upper left',
    frameon=False
)

#plt.tight_layout(rect=[0, 0, 0.85, 1])  # reserve space for legend on the right

plot_path = f'{prefix}/classification_compatible_mutations_{threshold}.pdf'
plt.savefig(plot_path, dpi=300, bbox_inches='tight')
print("Plot saved to:", plot_path)

