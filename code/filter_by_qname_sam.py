#!/usr/bin/env python3

import argparse
import pysam
from pathlib import Path

def load_qnames(files):
    q = set()
    for f in files:
        with open(f, "r") as fh:
            for line in fh:
                name = line.strip()
                if name:  # skip empty
                    q.add(name)
    return q

def guess_mode(path, write=False):
    """
    Return pysam mode for SAM/BAM based on extension.
    SAM: r/w ; BAM: rb/wb, CRAM not supported
    """
    ext = Path(path).suffix.lower()
    if write:
        return "wb" if ext == ".bam" else "w"
    else:
        return "rb" if ext == ".bam" else "r"

def main():
    ap = argparse.ArgumentParser(
        description="Filter SAM/BAM by QNAME using pysam (keep or drop)."
    )
    ap.add_argument("-i", "--in", dest="in_path", required=True, help="Input SAM/BAM")
    ap.add_argument("-o", "--out", dest="out_path", required=True, help="Output SAM/BAM")
    ap.add_argument("--qnames", nargs="+", required=True,
                   help="One or more files with QNAMEs (one per line, no header).")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--keep", action="store_true", help="Keep reads whose QNAME is in the list(s).")
    g.add_argument("--drop", action="store_true", default=False,
                   help="Drop reads whose QNAME is in the list(s).")
    ap.add_argument("--verbose", action="store_true", help="Print simple counters.")
    args = ap.parse_args()

    qnames = load_qnames(args.qnames)
    if args.verbose:
        print(f"Loaded {len(qnames):,} QNAMEs from {len(args.qnames)} file(s).")

    in_mode = guess_mode(args.in_path, write=False)
    out_mode = guess_mode(args.out_path, write=True)

    total = kept = dropped = 0
    with pysam.AlignmentFile(args.in_path, in_mode) as fin, \
         pysam.AlignmentFile(args.out_path, out_mode, header=fin.header) as fout:

        for aln in fin.fetch(until_eof=True):  # include unmapped/unsorted
            total += 1
            present = (aln.query_name in qnames)
            write_it = (present if args.keep else not present)
            if write_it:
                fout.write(aln)
                kept += 1
            else:
                dropped += 1

    if args.verbose:
        action = "kept" if args.keep else "dropped"
        print(f"Total: {total:,} | Kept: {kept:,} | Dropped: {dropped:,} | Mode: {action}")

if __name__ == "__main__":
    main()
