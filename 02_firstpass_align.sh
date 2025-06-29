#!/bin/bash
#SBATCH --job-name=Minimap2   # Job name
#SBATCH --mail-type=ALL       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sofia.kudasheva@earlham.ac.uk  # Where to send mail
#SBATCH --ntasks=1            # Run a single task
#SBATCH --cpus-per-task=12    # Number of CPU cores per task
#SBATCH --array=1-18          # Array range for 18 barcode folders
#SBATCH -p ei-medium           # Queue
#SBATCH --mem=60G            # Memory per CPU
#SBATCH --time=2-00:45:00       # Adjust runtime as needed
#SBATCH -o /ei/projects/1/17339e77-1a5d-440a-abde-70e615a0f98b/data/logs/Minimap2_%A_%a.out # STDOUT
#SBATCH -e /ei/projects/1/17339e77-1a5d-440a-abde-70e615a0f98b/data/logs/Minimap2_%A_%a.err # STDERR

#####!!!!!Deactivate any conda envs running interactively before submitting to avoid library conflicts!!!!
source package 222eac79-310f-4d4b-8e1c-0cece4150333  #minimap2 2.24
source package 638df626-d658-40aa-80e5-14a275b7464b #samtools-1.15.1

runs=(1 2)
# Define the base directory for the barcode folders
base_dir="/ei/projects/1/17339e77-1a5d-440a-abde-70e615a0f98b/data/results/20241031_IPSCs/processed_fastq/_${runs[1]}"
output_dir="/ei/projects/1/17339e77-1a5d-440a-abde-70e615a0f98b/data/results/20241031_IPSCs/alignments/_${runs[1]}"
genome="/ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/references/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
ref_juncs="/ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/references/gencode.v46.annotation.annotated_juncs.bed"
twopasstools="/ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/tools/2passtools/2passtools.sif"

# Create an array of barcode folder names
barcode_list=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18)
barcode=${barcode_list[$SLURM_ARRAY_TASK_ID-1]}
input_fastq="${base_dir}/barcode${barcode}/merged.processed.restranded.fastq.gz"  # Adjust if filenames are different
output_bam="${output_dir}/barcode${barcode}.sorted.bam"
echo "Processing ${input_fastq}, writing to ${output_bam}"

cpus=$SLURM_CPUS_PER_TASK

# # Run minimap2
minimap2 -ax splice --cs=long -I 16G -G 500k -ub \
        -k 14 --secondary=yes -t $cpus $genome $input_fastq | samtools view -b | samtools sort -O BAM -m 3G -@ $cpus -o $output_bam
        samtools index $output_bam
echo "Finished processing ${barcode}"

singularity exec --env LC_ALL=C.UTF-8,LANG=C.UTF-8 $twopasstools 2passtools score -f $genome -p $cpus -a $ref_juncs -m "GTAG|GCAG|ATAG|ATAC" -j 3 \
        --keep-all-annot -o ${output_dir}/barcode${barcode}.juncs.bed $output_bam


