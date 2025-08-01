#!/bin/bash
#SBATCH --job-name=Minimap2   # Job name
#SBATCH --mail-type=ALL       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sofia.kudasheva@earlham.ac.uk  # Where to send mail
#SBATCH --ntasks=1            # Run a single task
#SBATCH --cpus-per-task=12    # Number of CPU cores per task
#SBATCH --array=1-18          # Array range for 18 barcode folders
#SBATCH -p ei-medium           # Queue
#SBATCH --mem=80G            # Memory per CPU
#SBATCH --time=1-00:45:00       # Adjust runtime as needed
#SBATCH -o /ei/projects/1/17339e77-1a5d-440a-abde-70e615a0f98b/data/logs/Minimap2_%A_%a.out # STDOUT
#SBATCH -e /ei/projects/1/17339e77-1a5d-440a-abde-70e615a0f98b/data/logs/Minimap2_%A_%a.err # STDERR

source package 222eac79-310f-4d4b-8e1c-0cece4150333  #minimap2 2.24
source package 638df626-d658-40aa-80e5-14a275b7464b #samtools-1.15.1

base_dir="/ei/projects/1/17339e77-1a5d-440a-abde-70e615a0f98b/data/results/20241031_IPSCs/processed_fastq/merged_barcodes"
output_dir="/ei/projects/1/17339e77-1a5d-440a-abde-70e615a0f98b/data/results/20241031_IPSCs/alignments/merged_barcodes"
genome="/ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/references/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
filt_juncs="/ei/projects/1/17339e77-1a5d-440a-abde-70e615a0f98b/data/results/20241031_IPSCs/alignments/merged_juncs.filt.bed"

# Create an array of barcode folder names
barcode_list=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18)
barcode=${barcode_list[$SLURM_ARRAY_TASK_ID-1]}
input_fastq="${base_dir}/barcode${barcode}.fastq.gz" 
output_bam="${output_dir}/barcode${barcode}.sorted.bam"

cpus=$SLURM_CPUS_PER_TASK

# # Run minimap2
minimap2 -ax splice -k 14 -t $cpus -ub -G 500k --junc-bed $filt_juncs --junc-bonus 15 --secondary=yes \
        $genome $input_fastq | samtools sort -O BAM -m 3G -@ $cpus -o $output_bam
samtools index $output_bam
echo "Finished processing ${barcode}"


