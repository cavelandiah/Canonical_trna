#!/usr/bin/env python

import pandas as pd
import numpy as np
import re
import sys

sample = sys.argv[1]
threshold = sys.argv[2]
path_file = f"./Results/concatenated_{sample}_{threshold}.csv"
# Read just the header first
all_cols = pd.read_csv(path_file, sep=",", nrows=0,dtype=str, low_memory=False).columns
# Pick the ones we want
usecols = ["Accessibility"] + [c for c in all_cols if c.endswith("_vec")]

# Now read only those
df = pd.read_csv(path_file, sep=",", usecols=usecols, low_memory=False)
df['Sample'] = sample

# ---------- 2) Identify vec columns ----------
vec_cols = [c for c in df.columns if re.search(r"_vec$", c)]
# ---------- 3) Masks ----------
acc = df["Accessibility"].astype(int)
allowed_mask = acc == 1            # allowed positions
non_allowed_mask = ~allowed_mask   # non-allowed positions

allowed_n = int(allowed_mask.sum())
non_allowed_n = int(non_allowed_mask.sum())

# ---------- 4) Summaries per vec ----------
rows = []
for col in vec_cols:
    vec = df[col].fillna(0).astype(int)

    # Confusion-like counts relative to Accessibility
    TP = int(((vec == 1) & allowed_mask).sum())        # correct 1 in allowed
    FN = int(((vec == 0) & allowed_mask).sum())        # miss: 0 in allowed
    FP = int(((vec == 1) & non_allowed_mask).sum())    # violation: 1 in non-allowed
    TN = int(((vec == 0) & non_allowed_mask).sum())    # correct 0 in non-allowed

    # Normalized rates (guard divide-by-zero)
    tpr_allowed = TP / allowed_n if allowed_n else np.nan   # sensitivity on allowed=1 set
    fpr_non_allowed = FP / non_allowed_n if non_allowed_n else np.nan

    # Grouping rule (two groups)
    group = "violates_non_allowed" if FP > 0 else "allowed_only"

    # Optional stricter group: perfect everywhere
    perfect = (FP == 0 and FN == 0)

    rows.append({
        "vec": col,
        "TP_allowed_1match": TP,
        "FN_allowed_0mismatch": FN,
        "FP_nonallowed_1violation": FP,
        "TN_nonallowed_0match": TN,
        "TPR_allowed": tpr_allowed,
        "FPR_non_allowed": fpr_non_allowed,
        "allowed_sites": allowed_n,
        "non_allowed_sites": non_allowed_n,
        "group": group,
        "perfect_match": perfect,
        "sample": sample
    })

summary_df = pd.DataFrame(rows).sort_values(["group","vec"]).reset_index(drop=True)
print(summary_df)
summary_df.to_csv(f'./Results/classification_reads_{sample}_{threshold}.csv', index=False)
