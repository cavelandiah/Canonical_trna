#!/usr/bin/env python3

import RNA
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import varnaapi
import random
import shutil

from varnaapi.param import BasesStyle

## varnaapi conf
path_to_VARNA=shutil.which("VARNAv3-93.jar")
varnaapi.set_VARNA(path_to_VARNA)

def varna_annotate_profile(v, profile,
    colormap = [
        (1.2, "#FF6655"),
        (0.7, "#CC3300"),
        (0.3, "#DDCC00"),
        (0.0, "#000000"),
        (-0.5, "#4444BB"),
        ]
    ):
    """
    Annotate colors to varna structure plot based on a reactivity profile. The default profile works for SHAPE.
    Note: the profile vector is 0-based
    """

    # collect positions for each color
    collists = [[] for _ in colormap]
    for i in range(0,len(sequence)):
        for (t,c),l in zip(colormap,collists):
            if profile[i]>=t:
                l.append(i+1)
                break

    # set the colors
    white = "#FFFFFF"
    black = "#000000"
    for (t,c),positions in zip(colormap,collists):
        style = BasesStyle(fill=c, outline=c,
                           label=black if c!=black else white)
        v.add_bases_style(style, positions)

    return v

def plot_varna_react(seq, reactivities, fc, name='example'):
    fc.sc_add_SHAPE_deigan(reactivities, m=1.6, b=-0.6)
    ss, mfe = fc.mfe()
    print(f"SHAPE: {ss}, {mfe}")
    v = varnaapi.Structure(structure=ss, sequence=seq)
    v = varna_annotate_profile(v, reactivities)
    v.savefig(f"{name}_shape.eps")
    #v.show()
    return

#>iMet-CAT-2-1(iMet-CAT)
#AGCAGAGTGGCGCAGCGGAAGCGTGCTGGGCCCATAACCCAGAGGTCGATGGATCTAAACCATCCTCTGCTACCA
sequence = "AGCAGAGTGGCGCAGCGGAAGCGTGCTGGGCCCATAACCCAGAGGTCGATGGATCTAAACCATCCTCTGCTACCA"
# Reactivities
react_file="./reac.txt"
reactivities = []
with open(react_file, 'r') as fn:
    lines = fn.readlines()
reactivities = [float(line.strip()) for line in lines]


### Model details
md = RNA.md()
md.uniq_ML = 1
#md.temperature = 37.0
###
# Create fold compound
fc = RNA.fold_compound(sequence, md)
(ss, mfe) = fc.mfe()
print(f"BASIC: {ss}, {mfe}")

# Estimate the scale of the partition function, from MFE.
# This is used to rescale pf computation.
fc.exp_params_rescale(mfe)
# SS plot
plot_varna_react(sequence, reactivities, fc, name='iMet-CAT-2-1')
