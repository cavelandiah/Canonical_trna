#!/usr/bin/env python3

import sys
import pandas as pd

# Defs
def clean_vec_to_int(s: pd.Series) -> pd.Series:
    """Clean values: strip spaces, drop suffix, convert to integer indices"""
    s = s.astype(str).str.strip().str.replace(r"_vec$", "", regex=True)
    s = pd.to_numeric(s, errors="coerce").dropna().astype(int)
    return s

# Input
sample = sys.argv[1]
threshold = sys.argv[2]
test = int(sys.argv[3])

if test == 1:
    input_folder = "../results/Test"
else:
    input_folder = "../results"

# Read inputs
summary_df = pd.read_csv(f'{input_folder}/classification_reads_{sample}_{threshold}.csv')
df_index = pd.read_csv(f"{input_folder}/index_{sample}_{threshold}.csv", sep=r"\s+", engine="python")

# Subset columns
required_summary_cols = {"group", "vec"}
required_index_cols   = {"seq_idx", "QNAME"}
missing_summary = required_summary_cols - set(summary_df.columns)
missing_index   = required_index_cols   - set(df_index.columns)
if missing_summary:
    raise ValueError(f"Missing columns in summary_df: {missing_summary}")
if missing_index:
    raise ValueError(f"Missing columns in df_index: {missing_index}")

# select vec by group
compatible_summary = summary_df.loc[summary_df["group"] == "Compatible", "vec"]
noncompatible_summary = summary_df.loc[summary_df["group"] == "Noncompatible", "vec"]
noninformative_summary = summary_df.loc[summary_df["group"] == "Non-informative", "vec"]

compat_idx = clean_vec_to_int(compatible_summary)
noncompat_idx = clean_vec_to_int(noncompatible_summary)
noninfcompat_idx = clean_vec_to_int(noninformative_summary)

df_index = df_index.copy()
df_index["seq_idx"] = pd.to_numeric(df_index["seq_idx"], errors="coerce").astype("Int64")

comp_names = df_index.loc[df_index["seq_idx"].isin(compat_idx), "QNAME"]
noncomp_names = df_index.loc[df_index["seq_idx"].isin(noncompat_idx), "QNAME"]
noninf_names = df_index.loc[df_index["seq_idx"].isin(noninfcompat_idx), "QNAME"]

df_index_compatible = f"{input_folder}/index_{sample}_{threshold}_compatible.csv"
df_index_noncompatible = f"{input_folder}/index_{sample}_{threshold}_noncompatible.csv"
df_index_noninformative = f"{input_folder}/index_{sample}_{threshold}_noninformative.csv"

comp_names.to_csv(df_index_compatible, index=False, header=False)
noncomp_names.to_csv(df_index_noncompatible, index=False, header=False)
noninf_names.to_csv(df_index_noninformative, index=False, header=False)
