#!/usr/bin/env python3

import pysam
import os
import sys

sample= sys.argv[1]
thres = sys.argv[2]
sam=os.path.join('Data',str(sample)+".sam")
outsam=os.path.join('Data',str(sample)+f"_{thres}.sam")

# Open the input SAM file
with pysam.AlignmentFile(sam, "r") as infile, pysam.AlignmentFile(outsam, "w", header=infile.header) as outfile:
    for read in infile:
        if not read.is_unmapped:  # skip unmapped reads
            read_length = read.query_length  # length of the read sequence
            if read_length > int(thres):
                outfile.write(read)
