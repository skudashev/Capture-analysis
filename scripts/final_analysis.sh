#!/bin/bash
#SBATCH --job-name=Minimap2   # Job name
#SBATCH --mail-type=ALL       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sofia.kudasheva@earlham.ac.uk  # Where to send mail
#SBATCH --ntasks=1            # Run a single task
#SBATCH --cpus-per-task=8    # Number of CPU cores per task
#SBATCH -p ei-short           # Queue
#SBATCH --mem=50G            # Memory per CPU
#SBATCH --time=12:00:00       # Adjust runtime as needed
#SBATCH -o /ei/projects/8/8289c66d-2d56-4706-a307-5a9a3eb3747e/data/Capture/merged_capture/FINAL_RESULTS_March_2025/logs/Minimap2_%A.out # STDOUT
#SBATCH -e /ei/projects/8/8289c66d-2d56-4706-a307-5a9a3eb3747e/data/Capture/merged_capture/FINAL_RESULTS_March_2025/logs/Minimap2_%A.err # STDERR

source package 222eac79-310f-4d4b-8e1c-0cece4150333  #minimap2 2.24
source ~/.bashrc
conda activate lr_analysis

# Define directories and resources
output_dir="/ei/projects/8/8289c66d-2d56-4706-a307-5a9a3eb3747e/data/Capture/merged_capture/FINAL_RESULTS_March_2025"
cpus=$SLURM_CPUS_PER_TASK

# Input files
input_fastq="/ei/projects/8/8289c66d-2d56-4706-a307-5a9a3eb3747e/data/Capture/merged_capture/merged_capture-restranded.fastq.gz"
# Output files
out_bam="${output_dir}/merged.secondpass.sorted.bam"
filtered_bam="${output_dir}/merged.secondpass.filt.bam"
out_txt="${output_dir}/reads_to_filter.txt"

# References
genome="/ei/projects/3/31655266-640a-41d2-8663-59bba38bc3c4/data/data/References/hg38_sequin.fa"
ref_juncs="/ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/references/gencode.v46.annotation.annotated_juncs.bed"
# Tools
twopasstools="/ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/tools/2passtools/2passtools.sif"
samtools="/ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/tools/samtools/samtools-1.20.img"
scripts="/ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/scripts"

# Step 1: Apply 2passtools filtering
echo "Running 2passtools filtering..."
# singularity exec --env LC_ALL=C.UTF-8,LANG=C.UTF-8 $twopasstools 2passtools filter --exprs \
#     'count > 3 and jad > 3 and ((acceptor_seq_score > 0.5 and donor_seq_score > 0.5) or decision_tree_1_pred or decision_tree_2_pred)' \
#     -o ${output_dir}/merged.juncs.filt.bed ${output_dir}/merged.juncs.all.bed

# Step 2: Align with minimap2 and filter through samtools directly
echo "Running minimap2 alignment and sorting..."
# minimap2 -ax splice -k 14 -t $cpus --junc-bed ${output_dir}/merged.juncs.filt.bed --junc-bonus 12 \
#     --secondary=yes -ub -G 500000 $genome $input_fastq | \
#     singularity exec $samtools samtools sort -O BAM -m 3G -@ $cpus -o $out_bam
# singularity exec $samtools samtools index -@ $cpus $out_bam

# # Step 3: Filter reads using Python script
# echo "Running read filtering..."
# python $scripts/filter_RT.py --in_bam $out_bam --out_txt $out_txt

# # Step 4: Check if filtering file exists and is not empty
# if [ ! -s $out_txt ]; then
#     echo "No reads to filter. Exiting."
#     exit 1
# else
#     echo "Applying read filtering..."
#     singularity exec $samtools samtools view -@ $cpus -bh -N ^$out_txt -F 0x804 $out_bam | \
#     singularity exec $samtools samtools sort -@ $cpus -o $filtered_bam

#     singularity exec $samtools samtools index -@ $cpus $filtered_bam
# fi


# echo "Splitting BAM into chromosome-specific SAM files..."

# # Generate .fai index for genome if not already present
# if [ ! -f "${genome}.fai" ]; then
#     echo "Indexing genome..."
#     singularity exec $samtools samtools faidx $genome
# fi

# # Create output directory for chromosome SAM files
# chr_sam_dir="${output_dir}/per_chromosome_sam"
# mkdir -p $chr_sam_dir

# # Loop through each chromosome in the reference
# cut -f1 ${genome}.fai | grep -v '_' | while read chr; do
#     echo "Extracting reads for chromosome: $chr"
#     singularity exec $samtools samtools view -@ $cpus -h $filtered_bam $chr > "${chr_sam_dir}/${chr}.sam"
# done

python $scripts/split_bam_by_windows.py merged.secondpass.filt.bam merged.secondpass.filt 10000000

