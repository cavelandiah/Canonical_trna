#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
tRNA secondary structure with/without SHAPE (Deigan) using ViennaRNA Python API,
and colored rendering via VARNA.

Usage:
  ./generate_ss_shape_V2.py <selection_origin> <subsetting> <tto>
    selection_origin: 'g' or 'm1g'
    subsetting: '' or any suffix (e.g., 'subsetA')
    tto: 'NAI' or 'DMS'

Example:
  ./fold_trna_shape.py g '' NAI
"""

import sys
import os
import math
import shutil
import numpy as np
import pandas as pd
import RNA            # ViennaRNA Python bindings
import varnaapi
from varnaapi.param import BasesStyle

# ---------------------------
# Configuration (edit below)
# ---------------------------

# Sequence under study (iMet-CAT-2-1 as in your example)
SEQUENCE = "AGCAGAGTGGCGCAGCGGAAGCGTGCTGGGCCCATAACCCAGAGGTCGATGGATCTAAACCATCCTCTGCTACCA"

# Reactivity files live here; one value per line; use -999 or NaN for missing.
REACT_FOLDER = "./Reactivities"

# Default Deigan parameters (commonly used for SHAPE-MaP)
DEIGAN_M = 1.6
DEIGAN_B = -0.6

# VARNA jar autodetect (or hardcode full path)
PATH_TO_VARNA = shutil.which("VARNAv3-93.jar")

# SHAPE color scale (thresholds descending)
# Use a monotone mapping; we’ll clip/round into these bins.
SHAPE_COLORS = [
    (1.2, "#FF6655"),  # high
    (0.7, "#CC3300"),
    (0.3, "#DDCC00"),
    (0.0, "#000000"),  # near-zero
    (-0.5, "#4444BB")  # special: map missing/negative sentinels here
]

# ---------------------------
# Utilities
# ---------------------------

def die(msg: str, code: int = 2):
    print(f"[error] {msg}", file=sys.stderr)
    sys.exit(code)

def ensure_varna():
    if PATH_TO_VARNA is None:
        die("Could not find VARNAv3-93.jar on PATH. Please install/point PATH_TO_VARNA.")
    varnaapi.set_VARNA(PATH_TO_VARNA)

def build_names(selection_origin: str, subsetting: str, tto: str, lenrna: int):
    selection_origin = selection_origin.lower()
    tto = tto.upper()
    if selection_origin not in {"g", "m1g"}:
        die("selection_origin must be 'g' or 'm1g'.")

    if tto not in {"NAI", "DMS"}:
        die("tto must be 'NAI' or 'DMS'.")

    if selection_origin == "g":
        control = "G_i1"
        comparison = "G_i2" if tto == "NAI" else "G_i3"
    else:  # m1g
        control = "G_i4"
        comparison = "G_i6" if tto == "NAI" else "G_i7"

    if subsetting:
        control = f"{control}_{lenrna}_{subsetting}"
        comparison = f"{comparison}_{lenrna}_{subsetting}"
    else:
        control = f"{control}_{lenrna}"
        comparison = f"{comparison}_{lenrna}"

    # Original code had an uppercase TTO in the path; keep consistent & safe
    react_dir = os.path.join(REACT_FOLDER, tto)
    react_file = os.path.join(
        react_dir,
        r"iMet-CAT-2-1(iMet-CAT)_" + f"{control}-{comparison}.varna"
    )
    return control, comparison, react_file

def load_reactivities(path: str, n: int):
    if not os.path.isfile(path):
        die(f"Reactivity file not found: {path}")

    vals = []
    with open(path, "r") as fh:
        for line in fh:
            s = line.strip()
            if s == "":
                vals.append(np.nan)
            else:
                try:
                    vals.append(float(s))
                except ValueError:
                    # tolerate junk: mark as missing
                    vals.append(np.nan)

    # Pad or trim to match sequence length
    if len(vals) < n:
        vals += [np.nan] * (n - len(vals))
    elif len(vals) > n:
        vals = vals[:n]

    # Replace NaN with sentinel used by ViennaRNA's SHAPE parsing (-999)
    vals = [(-999.0 if (v is np.nan or (isinstance(v, float) and math.isnan(v))) else float(v)) for v in vals]
    return vals

def apply_shape_deigan(fc: RNA.fold_compound, reactivities, m=DEIGAN_M, b=DEIGAN_B):
    """
    Attach Deigan pseudo-energies to the fold compound.
    Missing values must be -999; ViennaRNA ignores them.
    """
    # sc_add_SHAPE_deigan returns 0 on failure; raise if that happens
    ok = fc.sc_add_SHAPE_deigan(reactivities, m=m, b=b)
    if ok == 0:
        die("sc_add_SHAPE_deigan() failed. Check reactivities length/sentinels.")
    return fc

def mfe_and_pf(fc: RNA.fold_compound):
    """
    Compute MFE and then PF on the SAME parameterization (important).
    Rescale PF using MFE estimate to stabilize numerics for long sequences.
    """
    ss_mfe, e_mfe = fc.mfe()
    fc.exp_params_rescale(e_mfe)  # rescale before partition function
    e_pf = fc.pf()                # returns PF free energy
    # centroid requires PF to be computed
    ss_centroid, dist = fc.centroid()
    return (ss_mfe, e_mfe), (ss_centroid, dist), e_pf

def annotate_varna_by_profile(v, profile, colormap):
    """
    Color nucleotides according to profile thresholds.
    Indices for VARNA are 1-based.
    We also map -999 (missing) to the last color band (threshold -0.5).
    """
    # Prepare bins
    bins = [[] for _ in colormap]
    for i, val in enumerate(profile):
        x = val
        # Treat missing sentinel explicitly
        if val == -999.0:
            x = -0.5  # will fall into the last bin in our default map
        # Pick the first band whose threshold is <= x (since list is descending)
        for (thr, _), bucket in zip(colormap, bins):
            if x >= thr:
                bucket.append(i + 1)  # VARNA 1-based
                break

    white, black = "#FFFFFF", "#000000"
    for (thr, color), idx_list in zip(colormap, bins):
        style = BasesStyle(fill=color, outline=color,
                           label=(black if color != black else white))
        if idx_list:
            v.add_bases_style(style, idx_list)
    return v

def render_varna(seq: str, structure: str, profile, out_prefix: str):
    ensure_varna()
    v = varnaapi.Structure(structure=structure, sequence=seq)
    v = annotate_varna_by_profile(v, profile, SHAPE_COLORS)
    # Prefer vector output; pdf or svg are robust
    #v.savefig(f"{out_prefix}.pdf")
    # Optionally also export SVG for editing
    v.savefig(f"{out_prefix}.svg")

def print_summary(tag: str, mfe, centroid):
    (ss_mfe, e_mfe) = mfe
    (ss_centroid, mean_bp_dist) = centroid
    print(f"== {tag} ==")
    print(f"MFE       : {e_mfe:7.2f} kcal/mol")
    print(f"Structure : {ss_mfe}")
    print(f"Centroid  : {ss_centroid}  (mean bp distance = {mean_bp_dist:.2f})")
    #print(f"PF energy : {pf_energy:7.2f} kcal/mol")
    print()

# ---------------------------
# Main
# ---------------------------

def main():
    if len(sys.argv) != 5:
        die("Usage: ./fold_trna_shape.py <selection_origin: g|m1g> <subsetting: ''|suffix> <tto: NAI|DMS> <len:int>")

    selection_origin = sys.argv[1]
    subsetting = sys.argv[2]
    tto = sys.argv[3]
    lenrna = sys.argv[4]

    control, comparison, react_file = build_names(selection_origin, subsetting, tto, lenrna)
    seq = SEQUENCE
    n = len(seq)

    # Load SHAPE/NAI/DMS reactivities and sanitize to length n with -999 for missing
    reactivities = load_reactivities(react_file, n)

    # --- Baseline (no SHAPE) ---
    md = RNA.md()
    md.uniq_ML = True   # keep canonical ML handling (bool is clearer)
    md.dangles = 2      # default; explicit for reproducibility
    fc0 = RNA.fold_compound(seq, md)
    mfe0, centroid0, pf0 = mfe_and_pf(fc0)
    print_summary("BASELINE", mfe0, centroid0)

    # --- SHAPE-directed (Deigan) ---
    # Use a fresh fold_compound so the soft constraints don’t leak.
    fc = RNA.fold_compound(seq, md)
    apply_shape_deigan(fc, reactivities, m=DEIGAN_M, b=DEIGAN_B)
    mfe_s, centroid_s, pf_s = mfe_and_pf(fc)
    print_summary(f"SHAPE (Deigan m={DEIGAN_M}, b={DEIGAN_B})", mfe_s, centroid_s)

    # --- Render structures with VARNA (colored by SHAPE) ---
    out_base = f"iMet-CAT-2-1_{control}-{comparison}"
    render_varna(seq, mfe_s[0], reactivities, out_prefix=out_base + "_MFE_SHAPE")
    render_varna(seq, centroid_s[0], reactivities, out_prefix=out_base + "_CENTROID_SHAPE")

if __name__ == "__main__":
    main()

