#!/usr/bin/env python

import sys
import pandas as pd
import os
import pysam
import re
import numpy as np
from Bio import SeqIO

## Get canonical positions
# Based on the comparisons of the fasta sequences from tRNAs and the INFERNAL model,
# I retrieved the sequence K01390.1_442-514 to get canonical positions. At the same,
# the sequence X14835.1_6927-7002 is used to get the subscripted positions. At the end,
# we use both sequences to reconstruct the canonical reference, encoded in a dictionary,
# that makes the relation between positions in the alignment and corresponding canonical
# positions. This included the variable loop, and subscripted positions: 17a, 20a, and 20b.
# see a_to_c_dict.
# In the match_cm_position function, the alignment is curated, and each sequence append the
# 3'CCA, that is not included in the RFAM alignment.


# Data tables:
#G_i1.fastq.gz,iMetCATg,DMSO,0.5,iMet_G_DMSO,G_i1,atcgt,Human,1,Hsap,Control_Neg,NA
#----
#G_i2.fastq.gz,iMetCATg,NAI,0.5,iMet_G_NAI,G_i2,agcta,Human,1,Hsap,Control_Neg,NA
#G_i3.fastq.gz,iMetCATg,DMS,0.5,iMet_G_DMS,G_i3,gcata,Human,1,Hsap,Control_Neg,NA
#G_i4.fastq.gz,iMetCATm1g,DMSO,0.5,iMet_m1G_DMSO,G_i4,tctag,Human,1,Hsap,Control_Neg,NA
#----
#G_i6.fastq.gz,iMetCATm1g,NAI,0.5,iMet_m1G_NAI,G_i6,tagct,Human,1,Hsap,Control_Neg,NA
#G_i7.fastq.gz,iMetCATm1g,DMS,0.5,iMet_m1G_DMS,G_i7,actga,Human,1,Hsap,Control_Neg,NA


# Based on canonical structure, I defined a vector of candidate nucleotides
# that, in theory, should be exposed in a higher probability to the probing chemicals
# This is based on the non-paired positions on the canonical SS.
# THIS IS WITHOUT MODIFICATIONS
shape_nts_canonical = [0,8,9,14,15,16,17,"17a",18,19,20,"20a","20b",21,26,32,33,34,35,36,37,38,44,45,46,47,48,54,55,56,57,58,59,60,73,74,75,76]
# Include the modification at position 9
#shape_nts_canonical_imet = [0,8,9,14,15,16,17,"17a",18,19,20,"20a","20b",21,26,32,33,34,35,36,37,38,44,45,46,47,48,54,55,56,57,58,59,60,73,74,75,76]
## Convert to strings:
shape_nts_canonical = list(map(str, shape_nts_canonical))
#shape_nts_canonical_imet = list(map(str, shape_nts_canonical_imet))

# Automatically detect vector columns by excluding known non-vector columns
non_vector_columns = ['ref_pos', 'ref_nt', 'align_pos', 'canonical_pos', 'Accessibility', 'Experiment']

## Definitions

def load_sam_to_dataframe(sam_path, threshold):
    samfile = pysam.AlignmentFile(sam_path, "r")
    records = []
    for read in samfile.fetch():
        if read.query_length < int(threshold):
            continue
        records.append({
            "QNAME": read.query_name,
            "FLAG": read.flag,
            "RNAME": read.reference_name,
            "POS": read.reference_start + 1,
            "MAPQ": read.mapping_quality,
            "CIGAR": read.cigarstring,
            "SEQ": read.query_sequence
        })
    samfile.close()
    return pd.DataFrame(records)

def get_simple_pad(matched_seq: str, pos1_based: int,) -> str:
    return (" " * (pos1_based - 1)) + matched_seq

def get_padded_seq(matched_seq, pos, lead, tail, seq):
    """
    Add '>' '<' to characters and compare the matched string with the original sequence
    This defines a 'padded' sequence.
    """
    lead_seq = seq[0:tail].lower()
    tail_seq = seq[len(seq) - tail + 1:].lower()
    # insertions = len(seq)-lead-tail-len(matched_seq)

    padding_len = pos - 1

    diff = padding_len - len(lead_seq)

    if diff <= 0:
        lead_seq = lead_seq[(diff * (-1)):]
    else:
        lead_seq = " " * (diff - 1) + ">" + lead_seq

    padded_seq = lead_seq + matched_seq + tail_seq + "<"
    return padded_seq

def get_matched_seq(cigar, seq):
    """
    Based on cigar and seq from the read, get matching part of sequence and its length.
    Change mismatches to lowercase and keep matches to uppercase.
    """
    cigar = correct_cigar(cigar)
    # Based on cigar, get letters (block types)
    block_types = [i for i in cigar if not i.isdigit()]
    block_sizes = re.split("I|D|M|S|=|X", cigar)[0:-1]
    block_sizes = [int(s) for s in block_sizes]
    match_seq = ""
    l = 0  # length of ungapped matched seq
    for i, t in enumerate(block_types):
        if t == "S":
            seq = seq[block_sizes[i]:]
        elif t == "M":
            match_seq += seq[0: block_sizes[i]].upper()
            seq = seq[block_sizes[i]:]
            l += block_sizes[i]
        elif t == "=":
            match_seq += seq[0: block_sizes[i]].upper()
            seq = seq[block_sizes[i]:]
            l += block_sizes[i]
        elif t == "X":
            match_seq += seq[0: block_sizes[i]].lower()
            seq = seq[block_sizes[i]:]
            l += block_sizes[i]
        elif t == "I":
            seq = seq[block_sizes[i]:]
            #l += block_sizes[i]
        elif t == "D":
            match_seq += "-" * block_sizes[i]
    return match_seq, l

def correct_cigar(cigar):
    """
    context: It assumes that tRNA reference is composed as: tRNA+(N's). When added tail N's, segemehl
    reports cigar strings with mismatches (X) and not insertions (I) at the end of the sequence. Just
    replace last X by I.
    input: CIGAR string
    output: CIGAR modified string
    """
    if cigar.endswith('X'):
        new_cigar = cigar
        new_cigar = re.sub(r'X$', 'I', new_cigar)
        return new_cigar
    else:
        return cigar

def match_cm_position(alignment_fa):
    """
    Build dict with list of nt, len(seq), seq without '.'
    """
    handle = alignment_fa
    pos_dict = {}
    for record in SeqIO.parse(handle, "fasta"):
        rname = record.id.replace('/', '_')
        # Add here CCA, because ref seq contained CCA
        seq = str(record.seq) + "CCA"
        #seq = str(record.seq)
        # Enumerate those that are not matching '.' or '-' in the alignment vs. CM
        pos_list = [i for i,s in enumerate(seq) if s not in  '.-' ]
        pos_dict[rname] = (pos_list, len(seq), seq.replace('.', '').replace('-', ''))
        #pos_dict[rname] = (pos_list, len(seq), seq.replace('.', '').replace('_', ''))
    return pos_dict

def alignment_to_canonical_positions_mapper(alignment_fa):
    a_to_c_dict= {}

    #get subscripted positions
    #X14835.1/6927-7002 contains subscripted canonical positions at given zero based index
    pos_ref = match_cm_position(alignment_fa)['X14835.1_6927-7002'][0]
    a_to_c_dict[pos_ref[17]] = '17a'
    a_to_c_dict[pos_ref[21]] = '20a'
    a_to_c_dict[pos_ref[22]] = '20b'

    #get regular canonical positions
    #K01390.1_442-514 contains all canonical non subscripted positions
    pos_ref = match_cm_position(alignment_fa)['K01390.1_442-514'][0]
    variable_loop_start = 0
    variable_loop_end = 0
    for i, j in enumerate(pos_ref):
        a_to_c_dict[j] = i+1
        if i+1 == 45:
            variable_loop_start = j+1
        if i+1 == 46:
            variable_loop_end = j

    # annotate variable loop
    for j in range(variable_loop_start, variable_loop_end):
        a_to_c_dict[j] = 'e'
    return a_to_c_dict

def classify_vector(row, vector_column):
    """
    This function compares each nucleotide in the specified `vector_column` list with the reference nucleotide
    found in the `ref_nt` column of the input row. It classifies each element as:

    - None: if the value is missing (NaN) or an empty string,
    - 0: if the value matches the reference nucleotide (`ref_nt`),
    - 1: if the value differs from the reference nucleotide.
    """
    classification = {}
    for col in vector_column:
        if pd.isna(row[col]) or (row[col] == ""): # Empty
            classification[col] = None
        elif row[col] == row['ref_nt']: # Is a match -> 0
            classification[col] = 0
        else: # Error in seq: 1, it means nt is accessible
            classification[col] = 1
    return classification

# Input files
sample = sys.argv[1]
thres = sys.argv[2]
test = int(sys.argv[3])
ref_fasta ="../data/iMet.fa"
alignment_fa ="../data/tRNAs_aln.fa"

# Output folders
out_folder = "../results/"

# Test
# G_i1= iMet_G_DMSO
# G_i4= iMet_m1G_DMSO

if test == 1:
    input_folder='../data/test'
else:
    # Real data
    input_folder='../data/raw'

sam = os.path.join(f'{input_folder}',str(sample)+f"_{thres}.sam")
df = load_sam_to_dataframe(sam, thres)
ref_dict = SeqIO.to_dict(SeqIO.parse(ref_fasta, "fasta"))
# Match positions with the CM model
alignment_pos_dict =  match_cm_position(alignment_fa)
# Based on previous references, get canonical positions for all tRNAs dict
canonical_pos_dict = alignment_to_canonical_positions_mapper(alignment_fa)

## Read sam file
df_clean = df[df['RNAME'] == "iMet-CAT-2-1(iMet-CAT)"]
total_reads = len(df)
if total_reads > 0:
    groups = df_clean.groupby('RNAME')
    data = []
    for name, g_df in groups:
        if name != "iMet-CAT-2-1(iMet-CAT)":
            continue
        ref_seq = str(ref_dict[name].seq)
        ref_len = len(ref_seq)

        # Create a column 'mapped_seq' with the seq inferred from CIGAR
        g_df["mapped seq"] = g_df.apply(
            lambda row: get_matched_seq(row["CIGAR"], row["SEQ"])[0], axis=1
        )
        g_df['padded mapped seq'] = g_df.apply(
            lambda row: get_simple_pad(
                row["mapped seq"],
                row["POS"],
                ),
        # g_df["padded mapped seq"] = g_df.apply(
            # lambda row: get_padded_seq(
            # row["mapped seq"],
            # row["POS"],
            # 0,
            # 0,
            # row["SEQ"],
        # ),
        axis=1,
        )
        g_df['Experiment'] = sample
        seqs = g_df["padded mapped seq"].to_list()
        # Include spaces in relation to the reference sequence
        # and the padded mapped sequence
        seqs = [s + " " * (ref_len+1 - len(s)) for s in seqs]
        seqs = [list(s) for s in seqs]

        # List of list of nt from sequences
        alignment_df = pd.DataFrame(seqs)

        # Get read index
        mapping = pd.DataFrame({
            "seq_idx": alignment_df.index,
            "QNAME": g_df["QNAME"].values
            })
        if test == 0:
            mapping.to_csv(f"{out_folder}/index_{sample}_{thres}.csv", sep="\t", index=False)
        else:
            mapping.to_csv(f"{out_folder}/Test/index_{sample}_{thres}.csv", sep="\t", index=False)

canonical_ref = alignment_df.transpose().copy()
canonical_ref.reset_index(inplace=True)
canonical_ref.rename(columns = {'index':'ref_pos'}, inplace = True)
canonical_ref['ref_nt'] = canonical_ref.apply(lambda row: ref_seq[int(row['ref_pos'])] if int(row['ref_pos']) < len(ref_seq) else 'N', axis = 1)

# Translate to canonical structure
if name in alignment_pos_dict.keys():
    alignment_pos_ref = alignment_pos_dict[name][0]
    # Add +1 last position
    alignment_pos_ref += [alignment_pos_ref[-1]+1]
    # Compare with the canonical reference, add '.' for additional positions.
    canonical_pos_ref = [str(canonical_pos_dict[pos])
                            if pos in canonical_pos_dict.keys()
                            else '.'
                            for pos in alignment_pos_ref]

# Add the column describing alignment position of the nt
canonical_ref['align_pos'] = canonical_ref.apply(lambda row: alignment_pos_ref[int(row['ref_pos'])], axis = 1)
# Add the column describing the canonical position
canonical_ref['canonical_pos'] = canonical_ref.apply(lambda row: canonical_pos_ref[int(row['ref_pos'])],axis = 1)
# Calculate Accessibility: If the nucleotide is or not accessible (based on the canonical structure: 1 is accessible, 0 paired.
# Here is a problem when the canonical position is like '.', I didn't enconded that.
# Take into account this:
# >iMet-CAT-2-1(iMet-CAT)
#AGCAGAGU..G..GCGCAGC...GGA...AGCGUG..C.UGGGC.CCAUAACCCA.GAG.
#..........................G.uCGAU.GGA.UCU.AA..ACCAU.C.CUC...
#UGCU.A
canonical_ref['Accessibility'] = canonical_ref['canonical_pos'].apply(lambda x: 1 if str(x) in shape_nts_canonical else 0)
canonical_ref.replace(" ", pd.NA, inplace=True)

# Automatically detect vector columns by excluding known non-vector columns
vector_columns = [col for col in canonical_ref.columns if col not in non_vector_columns]
# Apply the vector classification to each row
classified_vectors = canonical_ref.apply(lambda row: classify_vector(row, vector_columns), axis=1)
# Convert the resulting series of dictionaries into a DataFrame
classified_df = pd.DataFrame(list(classified_vectors))
classified_df = classified_df.add_suffix('_vec')
# Merge classified results back into the original DataFrame and remove last row
canonical_table_classified = pd.concat([canonical_ref, classified_df], axis=1).iloc[:-1]
# Write to Results:
if test == 0:
    canonical_table_classified.to_csv(f'{out_folder}/concatenated_{sample}_{thres}.csv',
                                      na_rep="NA",index=False)
else:
    canonical_table_classified.to_csv(f'{out_folder}/Test/concatenated_{sample}_{thres}.csv',
                                      na_rep="NA",index=False)
