#!/usr/bin/env python3
import os
import sys

if len(sys.argv) != 3:
    print("Usage: python make_tama_filelist.py <collapse_dir> <output_file>")
    sys.exit(1)

collapse_dir = sys.argv[1]
output_file = sys.argv[2]

with open(output_file, 'w') as out_f:
    for file in sorted(os.listdir(collapse_dir)):
        if file.endswith("collapsed.bed"):
            bed_path = os.path.join(collapse_dir, file)
            out_f.write(f"{bed_path}\tno_cap\t1,1,1\n")
