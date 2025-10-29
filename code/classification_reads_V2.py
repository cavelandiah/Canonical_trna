#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys

sample = sys.argv[1]
threshold = sys.argv[2]
test = int(sys.argv[3])

input_folder = "../results/"

if test == 1:
    path_file = f"{input_folder}/Test/concatenated_{sample}_{threshold}.csv"
else:
    # Real data
    path_file = f"{input_folder}/concatenated_{sample}_{threshold}.csv"

# --- 1) Read only required columns in one go ---
usecols = lambda c: c == "Accessibility" or c.endswith("_vec")
df = pd.read_csv(path_file, usecols=usecols, dtype=np.float32, low_memory=False)
df["Sample"] = sample

# --- 2) Pre-compute masks once ---
acc = pd.to_numeric(df["Accessibility"], errors="coerce").fillna(0).astype(np.int8)
allowed_mask = acc == 1
non_allowed_mask = ~allowed_mask

allowed_n = int(allowed_mask.sum())
non_allowed_n = len(acc) - allowed_n  # faster than non_allowed_mask.sum()

# --- 3) Collect vector columns ---
vec_cols = [c for c in df.columns if c.endswith("_vec")]

# --- 4) Pre-convert all vec columns to int8 matrix ---
X = df[vec_cols].apply(pd.to_numeric, errors="coerce").fillna(0).astype(np.int8)

# --- 5) Compute confusion-like metrics vectorized ---
TP = (X.T @ allowed_mask.values).astype(int)
FN = (X.shape[0] - X.T @ allowed_mask.values).astype(int)  # we'll fix this below
# Correction: FN counts 0s in allowed → sum((X==0) & allowed_mask)
FN = ((1 - X).T @ allowed_mask.values).astype(int)
FP = (X.T @ non_allowed_mask.values).astype(int)
TN = ((1 - X).T @ non_allowed_mask.values).astype(int)

# --- 6) Normalize and assemble results ---
summary_df = pd.DataFrame({
    "vec": vec_cols,
    "TP_allowed_1match": TP,
    "FN_allowed_0mismatch": FN,
    "FP_nonallowed_1violation": FP,
    "TN_nonallowed_0match": TN,
    "TPR_allowed": TP / allowed_n if allowed_n else np.nan,
    "FPR_non_allowed": FP / non_allowed_n if non_allowed_n else np.nan,
    "allowed_sites": allowed_n,
    "non_allowed_sites": non_allowed_n,
    "group": np.where(FP > 0, "violates_non_allowed", "allowed_only"),
    "perfect_match": (FP == 0) & (FN == 0),
    "sample": sample
})

summary_df = summary_df.sort_values(["group", "vec"]).reset_index(drop=True)
if test == 1:
    summary_df.to_csv(f"{input_folder}/Test/classification_reads_{sample}_{threshold}.csv", index=False)
else:
    summary_df.to_csv(f"{input_folder}/classification_reads_{sample}_{threshold}.csv", index=False)
