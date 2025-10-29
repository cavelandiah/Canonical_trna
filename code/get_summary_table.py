#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from scipy.spatial.distance import hamming, cosine


## Defs
# def n_lower_chars(string):
    # return sum(1 for c in string if c.islower())

# def n_ones_count(s):
    # """Count number of '1' = mutations in read
    # """
    # # Numpy array
    # return int(np.sum(s))

# def hamming_similarity(a, b):
    # """
    # Compute Hamming similarity between two binary vectors.
    # Counts matches at both 0 and 1 positions.

    # """
    # a = np.asarray(a)
    # b = np.asarray(b)
    # assert a.shape == b.shape, "Vectors must have the same length"

    # # Convert to strict binary in case of float values like 0.999999
    # a_bin = (a >= 0.5).astype(int)
    # b_bin = (b >= 0.5).astype(int)
    # return float(np.mean(a_bin == b_bin))

# def hamming_similarity_free_nt(ref, a):
    # """
    # Compute Hamming-like similarity restricted to positions where ref == 1.

    # Parameters
    # ----------
    # a : array-like
        # Read vector (0/1).
    # ref : array-like
        # Reference vector (0/1).

    # Returns
    # -------
    # float
        # Proportion of 1s in ref that are also 1 in a.
        # NaN if ref has no 1s.
    # """
    # a = np.asarray(a)
    # ref = np.asarray(ref)
    # assert a.shape == ref.shape, "Vectors must have same length"

    # mask_ref1 = (ref == 1)
    # if mask_ref1.sum() == 0:
        # return np.nan   # no 1s in ref to compare against

    # matches = np.sum(a[mask_ref1] == 1)
    # return matches / mask_ref1.sum()

# def hamming_distance_calc(ref, a):
    # a = np.asarray(a)
    # ref = np.asarray(ref)
    # hdistance = hamming(ref,a) #It returns normalized value
    # return hdistance


sample = sys.argv[1]
thres = sys.argv[2]
test = int(sys.argv[3])
if test == 0:
    canonical_table_classified = pd.read_csv(f'Results/concatenated_{sample}_{thres}.csv',
                                             low_memory=False)
    mapping = pd.read_csv(f"Results/index_{sample}_{thres}.csv", sep="\t")
else:
    mutation=sys.argv[4]
    canonical_table_classified = pd.read_csv(f'Data/Simulations/SAM/concatenated_{sample}_{thres}_{mutation}.csv',
                                             low_memory=False)
    mapping = pd.read_csv(f"Data/Simulations/SAM/index_{sample}_{thres}_{mutation}.csv", sep="\t")
    #canonical_table_classified = pd.read_csv(f'Results/Test/concatenated_{sample}_{thres}.csv',
    #                                         low_memory=False)
    #mapping = pd.read_csv(f"Results/Test/index_{sample}_{thres}.csv", sep="\t")


### Calculate distances table
# Extract vector columns once
vec_cols = [c for c in canonical_table_classified.columns if str(c).endswith("_vec")]
# Canonical vector
accessibility_vector = canonical_table_classified['Accessibility'].to_numpy(dtype=float)

read_ids, cos_dists, hamm_dists, free_precisions, n_mut = [], [], [], [], []

# map index->QNAME
idx2qname = dict(zip(mapping["seq_idx"].astype(str) + "_vec", mapping["QNAME"]))

for col in vec_cols:
    read_vec = canonical_table_classified[col].to_numpy()
    covered_mask = ~pd.isna(read_vec)
    if covered_mask.sum() == 0:
        continue

    read_sub = read_vec[covered_mask].astype(float, copy=False)
    acc_sub  = accessibility_vector[covered_mask].astype(float, copy=False)

    # binary projections
    read_bin = (read_sub >= 0.5).astype(int, copy=False)
    acc_bin  = (acc_sub  >= 0.5).astype(int, copy=False)

    # cosine distance
    #sim = cosine(acc_sub, read_bin).item()
    #if sim == 0:
    #    sim = np.nan
    #    cos_dists.append(np.nan)
    #else:
    #    cos_dists.append(sim)

    # Hamming distance (fraction disagreeing)
    #hamm_dists.append(float(hamming(acc_bin, read_bin)))

    # “free-nt” precision: among ref==1 positions, how many are 1 in read?
    m = (acc_bin == 1)
    free_precisions.append(np.nan if m.sum() == 0 else (read_bin[m].sum() / m.sum()))

    # mutation count across covered positions
    n_mut.append(int(read_bin.sum()))

    # read id
    read_ids.append(idx2qname.get(col, col))

distance_df = pd.DataFrame({
    'read': read_ids,
    #'cosine_distance': cos_dists,
    #'hamming_distance': hamm_dists,
    'precision_on_ref1': free_precisions,
    'mutations': n_mut,
    'Experiment': sample,
})

if test == 0:
    distance_df.to_csv(f'Results/distance_{sample}_{thres}.csv', na_rep='NA',index=False)
else:
    distance_df.to_csv(f'Data/Simulations/SAM/distance_{sample}_{thres}_{mutation}.csv', na_rep='NA',index=False)
    #distance_df.to_csv(f'Results/Test/distance_{sample}_{thres}.csv', na_rep='NA',index=False)

## OLD
# ## Calculare similarity table
# # Calculate similarity
# # Extract *_vec columns
# vec_cols = [str(col) for col in canonical_table_classified.columns if str(col).endswith('_vec')]
# # Convert Canonical Accessibility to binary vector
# accessibility_vector = canonical_table_classified['Accessibility'].to_numpy()

# read_ids = []
# similarities = []
# mutations = []
# hammingS = []
# hamming1 = []
# hammingD = []

# # Iteration over vector columns
# for col in vec_cols:
    # # Convert Canonical Accessibility to binary vector
    # read_vec = canonical_table_classified[col].to_numpy()
    # # Detects not NAs values and creates a vector with: True (not NA) or False (with NAs)
    # covered_mask = ~pd.isna(read_vec)
    # # Read should not have NAs, everything is True.
    # #if not covered_mask.all():
    # # Here check that at least 1 has coverage
    # if covered_mask.sum() == 0:
        # print(f"Read {col} has no coverage.")
        # continue

    # # With this, I select those that are True, it means select all that
    # # have coverage
    # read_subvec = read_vec[covered_mask].astype(float)
    # access_subvec = accessibility_vector[covered_mask].astype(float)

    # # Calculate distance
    # sim = cosine_similarity(read_subvec.reshape(1,-1),
                            # access_subvec.reshape(1, -1)).flatten()[0]
    # cos_d = 1 - sim
    # # Hamming
    # read_bin = (read_subvec >= 0.5).astype(int)
    # access_bin = (access_subvec >= 0.5).astype(int)
    # hamm = hamming_similarity(access_bin, read_bin)
    # hamm1 = hamming_similarity_free_nt(access_bin, read_bin)
    # hammD = 1 - hamm #hamming_distance_calc(access_bin, read_bin)
    # # Append results
    # read_ids.append(col)
    # similarities.append(cos_d)
    # hammingS.append(hamm)
    # hamming1.append(hamm1)
    # hammingD.append(hammD)
    # # Mutations per read
    # number_mutations = n_ones_count(read_subvec)
    # mutations.append(number_mutations)

# similarity_df = pd.DataFrame({
    # 'read': read_ids,
    # 'similarity_to_accessibility': similarities,
    # 'Hamming_similarity': hammingS,
    # 'Hamming_similarity_free': hamming1,
    # 'Hamming_distance': hammingD,
    # 'Mutations': mutations,
    # 'Experiment': sample,
    # })
# similarity_df.to_csv(f'Results/similarity_{sample}_{thres}.csv', na_rep='NA',index=False)
