#!/usr/bin/env python3
import sys
import os
import gffutils

classification_file = sys.argv[1]  # SQANTI classification tama_iq_transcriptome.filtered.probes.classification.txt
input_gtf = sys.argv[2]           # GTF  tama_iq_transcriptome.filtered.probes.gtf
reference_file = sys.argv[3]       # Reference GTF or .db /ei/projects/8/8289c66d-2d56-4706-a307-5a9a3eb3747e/data/Annotations/gencode.v40.annotation_sequin.gtf
output_file = sys.argv[4]          # Output file

# Only allow GTF input
if not input_gtf.endswith('.gtf'):
    print("Error: input file must be a GTF (.gtf)")
    sys.exit(1)

# Map SQANTI categories to suffixes
extension_dict = {
    "incomplete-splice_match": "ism",
    "novel_not_in_catalog": "nnc",
    "novel_in_catalog": "nic"
}

# Build replacement dictionary
replacement_dict = {}
with open(classification_file, 'r') as f:
    header = f.readline().strip().split('\t')
    for line in f:
        cols = line.strip().split('\t')
        transcript_id = cols[0]
        structural_category = cols[5]
        associated_transcript = cols[7]
        if structural_category in extension_dict:
            replacement_dict[transcript_id] = f"{transcript_id}.{extension_dict[structural_category]}"
        elif structural_category == "full-splice_match":
            replacement_dict[transcript_id] = f"{transcript_id}_{associated_transcript}"

# Prepare gffutils reference DB
if reference_file.endswith('.gtf'):
    db_file = reference_file.replace('.gtf', '.db')
elif reference_file.endswith('.db'):
    db_file = reference_file
else:
    print("Error: reference must be a .gtf or .db file")
    sys.exit(1)

if not os.path.exists(db_file) and reference_file.endswith('.gtf'):
    print("Creating gffutils database...")
    gffutils.create_db(reference_file, dbfn=db_file, force=True, keep_order=True,
                       merge_strategy='merge', sort_attribute_values=True,
                       disable_infer_transcripts=True, disable_infer_genes=True)

db = gffutils.FeatureDB(db_file, keep_order=True)

fsm_total = 0
fsm_matched = 0

# Process GTF input
with open(input_gtf, 'r') as f_in, open(output_file, 'w') as f_out:
    for line in f_in:
        if line.startswith('#') or not line.strip():
            f_out.write(line)
            continue

        fields = line.strip().split('\t')

        # Only process transcript features for coordinate matching
        if fields[2] != "transcript":
            f_out.write(line)
            continue

        # Extract transcript ID
        attrs = fields[8]
        tx_id = None
        for attr in attrs.split(';'):
            if 'transcript_id' in attr:
                tx_id = attr.split('"')[1]
                break

        if tx_id not in replacement_dict:
            f_out.write(line)
            continue

        new_tx_id = replacement_dict[tx_id]

        # FSM logic: check coords against reference
        if new_tx_id.endswith("_" + new_tx_id.split('_')[-1]):  # crude check for FSM
            fsm_total += 1
            try:
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]

                for feature in db.region(seqid=chrom, start=start, end=end, featuretype='transcript'):
                    if feature.strand == strand and feature.start == start and feature.end == end:
                        new_tx_id = f"{tx_id}_{feature.id}"
                        fsm_matched += 1
                        break
            except Exception as e:
                print(f"Warning: could not check coords for {tx_id}: {e}")

        # Replace transcript ID in the attributes field
        fields[8] = fields[8].replace(f'transcript_id "{tx_id}"', f'transcript_id "{new_tx_id}"')
        fields[8] = fields[8].replace(f'gene_id "{tx_id}"', f'gene_id "{new_tx_id}"')  # optional if gene_id matches tx_id

        f_out.write('\t'.join(fields) + '\n')

print(f"FSM transcripts total: {fsm_total}")
print(f"FSM transcripts matched reference coords: {fsm_matched}")
print(f"Output written to {output_file}")


#!/usr/bin/env python3
import pyranges as pr

def add_gene_ids_from_reference(input_gtf, ref_gtf, output_gtf):
    """
    Adds gene_id to input GTF based on gene_name mapping from reference GTF.

    Parameters:
        input_gtf (str): Path to the input GTF file.
        ref_gtf (str): Path to the reference GTF file.
        output_gtf (str): Path to the output GTF file.
    """
    
    # Read input and reference GTFs
    gr = pr.read_gtf(input_gtf)
    ref = pr.read_gtf(ref_gtf, full=True)
    
    # Extract reference genes
    ref_genes = ref.df[ref.df["Feature"] == "gene"].drop(ref.df.columns[16:], axis=1, errors="ignore")
    # Create mapping dict: gene_name → gene_id
    ref_name_to_id = {row["gene_name"]: row["gene_id"] for _, row in ref_genes.iterrows()}
    # Convert input PyRanges to DataFrame
    gr_df = gr.df.rename(columns={"gene_id": "gene_name"})
    
    # Add gene_id column based on mapping
    gr_df["gene_id"] = gr_df["gene_name"].apply(lambda x: ref_name_to_id.get(x))
    
    # Save back to PyRanges and write to GTF
    updated_gr = pr.PyRanges(gr_df).sort()
    updated_gr.to_gtf(output_gtf)
    
    print(f"Gene IDs added from reference and saved to {output_gtf}")

# Example usage
if __name__ == "__main__":
    import sys
    if len(sys.argv) != 4:
        print("Usage: python add_gene_ids.py <input.gtf> <reference.gtf> <output.gtf>")
        sys.exit(1)
    
    input_gtf = sys.argv[1]
    ref_gtf = sys.argv[2]
    output_gtf = sys.argv[3]
    
    add_gene_ids_from_reference(input_gtf, ref_gtf, output_gtf)
