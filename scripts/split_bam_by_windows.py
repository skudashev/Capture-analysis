import pysam
import sys
import os
import math

# Usage: python split_bam_by_windows.py input.bam output_prefix window_size

input_bam = sys.argv[1]
output_prefix = sys.argv[2]
window_size = int(sys.argv[3])  # e.g., 10000000 for 10 Mb

bam = pysam.AlignmentFile(input_bam, "rb")  # read BAM
header = bam.header

# Create output directory
output_dir = output_prefix + "_split"
os.makedirs(output_dir, exist_ok=True)

print(f"Splitting {input_bam} into {window_size:,} bp bins...")

for ref in bam.references:
    ref_len = bam.get_reference_length(ref)

    # skip non-standard chromosomes (with underscores)
    if "_" in ref:
        print(f"Skipping {ref}")
        continue

    num_bins = math.ceil(ref_len / window_size)

    for bin_index in range(num_bins):
        start = bin_index * window_size
        end = min((bin_index + 1) * window_size, ref_len)

        out_name = f"{output_dir}/{output_prefix}_{ref}_{start}_{end}.sam"
        out_sam = pysam.AlignmentFile(out_name, "w", header=header)

        print(f"Writing: {out_name}")
        for read in bam.fetch(ref, start, end):
            out_sam.write(read)

        out_sam.close()

bam.close()
print("Done.")
