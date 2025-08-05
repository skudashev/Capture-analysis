import pandas as pd
import argparse

# Parse the command line arguments
parser = argparse.ArgumentParser()
parser.add_argument("--sample_sheet", type=str, required=True, help="Path to the sample sheet file")
parser.add_argument("--input_dir", type=str, required=True, help="Path to the directory containing input BAM files")
parser.add_argument("--output_labels", type=str, required=True, help="Path to save the labels output")
args = parser.parse_args()

# Read the sample sheet
sample_sheet = pd.read_csv(args.sample_sheet, sep="\t", header=0, dtype=str)

# Generate BAM file paths based on sample_id
sample_sheet["BAM_File"] = sample_sheet["sample_id"].apply(
    lambda sid: f"{args.input_dir}/{sid}.filt.bam"
)

# Join all sample_id labels
labels = " ".join(sample_sheet["sample_id"])

# Save the labels to a file
with open(args.output_labels, "w") as f:
    f.write(labels + "\n")


