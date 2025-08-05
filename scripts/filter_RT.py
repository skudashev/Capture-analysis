import argparse
import pysam
import re

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Process a BAM file to find primary alignments with supplementary reads.")
    parser.add_argument('--in_bam', required=True, help="Input BAM file.")
    parser.add_argument('--out_txt', required=True, help="Output text file to save the filtered read IDs.")
    
    return parser.parse_args()

def process_bam(in_bam, out_txt):
    bam = pysam.AlignmentFile(in_bam, "rb")

    with open(out_txt, "w") as output_file:
        print("Processing BAM file for primary alignments with supplementary reads")

        # Process each read
        for read in bam.fetch():
            if read.has_tag('SA') and read.flag == 0:
                # Primary alignment (flag 0) with supplementary alignment
                ID = read.query_name
                query_s = read.query_alignment_start #index of the first base in query_sequence that is not soft-clipped.
                query_e = read.query_alignment_end 
                primary_length = query_e - query_s # Length of the primary alignment
                total_read_length = read.infer_read_length()

                # Check if primary alignment is less than 30% of the read length
                if primary_length < 0.3 * total_read_length:
                    print(ID, file=output_file)
                    print(ID + " is less than 30% of the read length")
                    continue  # No need to check overlap if this condition is met

                # Check for large deletion (D) over 20bp in the primary alignment's CIGAR string
                cigar_elements = read.cigartuples  # CIGAR string as a list of tuples (operation, length)
                for operation, length in cigar_elements:
                    # 2 represents a 'D' (deletion) in the CIGAR string
                    if operation == 2 and length > 20:
                        print(ID, file=output_file)
                        print(ID + " has a large deletion (>20bp) in the primary alignment")
                        continue # No need to check overlap if this condition is met

                primary_strand = '-' if read.is_reverse else '+'
                chr = read.reference_name
                ref_start = read.reference_start
                ref_end = read.reference_end
                
                # Check if supplementary alignment (SA tag) maps to reverse strand and overlaps
                supplementary_alignments = read.get_tag('SA').split(';')[:-1]
                for sa in supplementary_alignments:
                    sa_chr, sa_pos, sa_strand, sa_cigar, sa_mapq, sa_NM = sa.split(',')
                    if sa_chr == chr and sa_strand != primary_strand:
                        # Infer query length from CIGAR string of supplementary alignment
                        cigar_elements = re.findall(r'(\d+)([MDNX=])', sa_cigar)
                        sa_aligned_length = sum(int(length) for length, operation in cigar_elements if operation in 'MDNX=')
                        sa_pos = int(sa_pos)
                        if ref_start <= sa_pos <= ref_end or sa_pos <= ref_start <= sa_pos + sa_aligned_length:
                            print(ID, file=output_file)
                            print(ID + " has a supplementary alignment that overlaps with the primary alignment")
                            continue

    bam.close()

if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_arguments()

    # Process the BAM file with the provided arguments
    process_bam(args.in_bam, args.out_txt)


