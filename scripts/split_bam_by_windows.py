#!/usr/bin/env python3
import pysam
import os
import sys

if len(sys.argv) != 4:
    print("Usage: python split_bam_by_windows.py <input.bam> <window_size> <output_dir>")
    sys.exit(1)

input_bam = sys.argv[1]
window_size = int(sys.argv[2])
output_dir = sys.argv[3]

os.makedirs(output_dir, exist_ok=True)
bam = pysam.AlignmentFile(input_bam, "rb")

for ref in bam.header.references:
    length = bam.header.get_reference_length(ref)
    for start in range(0, length, window_size):
        end = min(start + window_size, length)
        out_bam_path = os.path.join(output_dir, f"{ref}_{start}_{end}.bam")
        with pysam.AlignmentFile(out_bam_path, "wb", header=bam.header) as out_bam:
            for read in bam.fetch(ref, start, end):
                out_bam.write(read)

bam.close()


