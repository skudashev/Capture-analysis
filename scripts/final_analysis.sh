#!/bin/bash
#SBATCH --job-name=Minimap2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@example.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH -p ei-short
#SBATCH --mem=50G
#SBATCH --time=12:00:00
#SBATCH -o logs/Minimap2_%A.out
#SBATCH -e logs/Minimap2_%A.err

# ================================
#        CONFIGURATION
# ================================

# Environment setup
source ~/.bashrc
conda activate lr_analysis

# Software packages (adjust as needed)
source package 222eac79-310f-4d4b-8e1c-0cece4150333   # minimap2 2.24
twopasstools_img="/path/to/2passtools.sif"
samtools_img="/path/to/samtools-1.20.img"

# Project directories
output_dir="/path/to/output"
mkdir -p "${output_dir}/logs" "${output_dir}/per_chromosome_sam"

# Inputs
input_fastq="/path/to/reads.fastq.gz"
genome="/path/to/genome.fa"
ref_juncs="/path/to/gencode_junctions.bed"

# Outputs
out_bam="${output_dir}/merged.secondpass.sorted.bam"
filtered_bam="${output_dir}/merged.secondpass.filt.bam"
out_txt="${output_dir}/reads_to_filter.txt"
juncs_filt_bed="${output_dir}/merged.juncs.filt.bed"

# Helper scripts
scripts_dir="/path/to/scripts"

cpus=$SLURM_CPUS_PER_TASK

# ================================
#        STEP 1: FILTER JUNCTIONS
# ================================
echo "[INFO] Running 2passtools filtering..."
singularity exec --env LC_ALL=C.UTF-8,LANG=C.UTF-8 "$twopasstools_img" \
  2passtools filter --exprs \
  'count > 3 and jad > 3 and ((acceptor_seq_score > 0.5 and donor_seq_score > 0.5) or decision_tree_1_pred or decision_tree_2_pred)' \
  -o "$juncs_filt_bed" "$output_dir/merged.juncs.all.bed"

# ================================
#        STEP 2: ALIGNMENT
# ================================
echo "[INFO] Running minimap2 alignment..."
minimap2 -ax splice -k 14 -t "$cpus" \
  --junc-bed "$juncs_filt_bed" --junc-bonus 12 \
  --secondary=yes -ub -G 500000 "$genome" "$input_fastq" | \
  singularity exec "$samtools_img" samtools sort -O BAM -m 3G -@ "$cpus" -o "$out_bam"

singularity exec "$samtools_img" samtools index -@ "$cpus" "$out_bam"

# ================================
#        STEP 3: FILTER READS
# ================================
echo "[INFO] Running read filtering..."
python "$scripts_dir/filter_RT.py" --in_bam "$out_bam" --out_txt "$out_txt"

if [ ! -s "$out_txt" ]; then
    echo "[INFO] No reads to filter. Skipping filtering step."
    cp "$out_bam" "$filtered_bam"
    singularity exec "$samtools_img" samtools index -@ "$cpus" "$filtered_bam"
else
    echo "[INFO] Filtering reads..."
    singularity exec "$samtools_img" samtools view -@ "$cpus" -bh -N ^"$out_txt" -F 0x804 "$out_bam" | \
      singularity exec "$samtools_img" samtools sort -@ "$cpus" -o "$filtered_bam"
    singularity exec "$samtools_img" samtools index -@ "$cpus" "$filtered_bam"
fi

# ================================
#        STEP 4: SPLIT BAM BY CHR
# ================================
echo "[INFO] Splitting BAM into chromosome-specific SAM files..."

if [ ! -f "${genome}.fai" ]; then
    echo "[INFO] Indexing genome..."
    singularity exec "$samtools_img" samtools faidx "$genome"
fi

chr_sam_dir="${output_dir}/per_chromosome_sam"
cut -f1 "${genome}.fai" | grep -v '_' | while read chr; do
    echo "[INFO] Extracting reads for chromosome: $chr"
    singularity exec "$samtools_img" samtools view -@ "$cpus" -h "$filtered_bam" "$chr" > "${chr_sam_dir}/${chr}.sam"
done

# ================================
#        STEP 5: SPLIT BAM BY WINDOWS
# ================================
echo "[INFO] Splitting BAM by windows..."
python "$scripts_dir/split_bam_by_windows.py" "$filtered_bam" "${output_dir}/split_bam" 10000000

echo "[INFO] Pipeline completed successfully."
