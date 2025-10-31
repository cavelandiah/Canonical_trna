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

# Get Accessibility and _vec columns
usecols = lambda c: c == "Accessibility" or c.endswith("_vec")
df = pd.read_csv(path_file, usecols=usecols, dtype=np.float32, low_memory=False)
df["Sample"] = sample

# df['Accessibility'] is float32 -> np.int8 (0,1)
acc = pd.to_numeric(df["Accessibility"], errors="coerce").fillna(0).astype(np.int8)
allowed_mask = acc == 1
non_allowed_mask = ~allowed_mask
allowed_n = int(allowed_mask.sum())
non_allowed_n = len(acc) - allowed_n

# Collect vector columns ---
vec_cols = [c for c in df.columns if c.endswith("_vec")]

# Pre-convert all vec columns to int8 matrix: on positions that are not numeric, or na,
# it converts to 0, in this analyses I do care about the positions of 1, other missing or 0s
# positions doesn't matter.
X = df[vec_cols].apply(pd.to_numeric, errors="coerce").fillna(0).astype(np.int8)

# Compute confusion-like metrics vectorized:
# TP: Number of positions that are free (1) in both: the canonical structure (allowed_mask) and in the read vector (X).
# FN: Number of positions that are free (1) in the canonical structure but not free (0) in the read vector (missed free positions).
# FP: Number of positions that aren't free (0) in the canonical structure (not_allowed_mask) and free in the read vector (X).
# TN: Number of positions that aren't free (0) in both: the canonical structure (not_allowed_mask) and free in the read vector (X).
# TPR_allowed measures how many canonical free positions were correctly identified as free in the vector.
# FPR_non_allowed measures how often the vector falsely marks a paired canonical position as free.
TP = (X.T @ allowed_mask.values).astype(int)
FN = ((1 - X).T @ allowed_mask.values).astype(int)
FP = (X.T @ non_allowed_mask.values).astype(int)
TN = ((1 - X).T @ non_allowed_mask.values).astype(int)

# Get summary
summary_df = pd.DataFrame({
    "vec": vec_cols,
    "TP": TP,
    "FN": FN,
    "FP": FP,
    "TN": TN,
    "TPR": TP / allowed_n if allowed_n else np.nan,
    "FPR": FP / non_allowed_n if non_allowed_n else np.nan,
    "allowed_sites": allowed_n,
    "non_allowed_sites": non_allowed_n,
    "group": np.where(FP > 0, "Noncompatible", "Compatible"),
    "perfect_match": (FP == 0) & (FN == 0),
    "sample": sample
})

summary_df = summary_df.sort_values(["group", "vec"]).reset_index(drop=True)
if test == 1:
    summary_df.to_csv(f"{input_folder}/Test/classification_reads_{sample}_{threshold}.csv", index=False)
else:
    summary_df.to_csv(f"{input_folder}/classification_reads_{sample}_{threshold}.csv", index=False)
