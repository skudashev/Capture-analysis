#!/usr/bin/env python3
import pandas as pd

# Input files
all_probes_file = "/ei/projects/8/8289c66d-2d56-4706-a307-5a9a3eb3747e/data/Capture/Capture_experiment_data/all_probes_names.txt"
lincRNA_file = "/ei/projects/8/8289c66d-2d56-4706-a307-5a9a3eb3747e/data/Capture/Capture_experiment_data/Capture_probes_list/lincRNA_probes.txt"
detected_genes_file = "detected_genes.txt"
output_file = "/ei/projects/8/8289c66d-2d56-4706-a307-5a9a3eb3747e/data/Capture/Capture_experiment_data/all_probes_annotated.txt"

# Load all probes (2 columns: symbol, gene_id)
df = pd.read_csv(all_probes_file, sep=r'\s+', header=None, names=["symbol", "gene_id"])

# Load lincRNA gene IDs
lincRNA_ids = pd.read_csv(lincRNA_file, header=None)[0].tolist()

# Load detected genes (symbols)
detected_symbols = pd.read_csv(detected_genes_file, header=None)[0].tolist()

# Function to classify category
def classify_category(symbol, gene_id):
    # Check by gene_id list
    if gene_id in lincRNA_ids:
        return "lnRNA"
    # Check by symbol patterns
    if symbol.startswith("ENSG"):
        return "lnRNA"
    if symbol.upper().startswith("LINC"):
        return "lnRNA"
    if "-AS" in symbol.upper() or "-DT" in symbol.upper():
        return "lnRNA"
    # Default
    return "protein_coding"

# Apply classification
df["category"] = df.apply(lambda row: classify_category(row["symbol"], row["gene_id"]), axis=1)

# Add detected/not_detected annotation
df["detection_status"] = df["symbol"].apply(lambda s: "detected" if s in detected_symbols else "not_detected")

# Save output
df.to_csv(output_file, sep="\t", index=False)

not_detected_counts = df[df["detection_status"] == "not_detected"]["category"].value_counts()
total_detected = (df["detection_status"] == "detected").sum()

# Print summary
print(f"Annotated file written to {output_file}")
print("\nSummary:")
print(f"Total detected genes: {total_detected}")
print(f"Protein coding not detected: {not_detected_counts.get('protein_coding', 0)}")
print(f"lnRNA not detected: {not_detected_counts.get('lnRNA', 0)}")


