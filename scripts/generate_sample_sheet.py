import pandas as pd

# Load your sample description file
df = pd.read_csv("Capture_sample_description.txt", sep="\t")

# Create the Path column based on Run and Barcode
df["Path"] = df.apply(
    lambda row: f"/ei/.project-scratch/4/461a2024-de28-4013-a1f0-549f4636fa28/Nicola_capture/Capture{row['Run']}/barcodes/barcode{str(row['Barcode']).zfill(2)}",
    axis=1
)

# rename Tube_ID to SampleID
df.rename(columns={"Tube_ID": "sample_id"}, inplace=True)

# Save as TSV
df.to_csv("samples.tsv", sep="\t", index=False)
