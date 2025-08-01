#!/bin/bash
#SBATCH --job-name=process_fq   # Job name
#SBATCH --mail-type=ALL       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sofia.kudasheva@earlham.ac.uk  # Where to send mail
#SBATCH --ntasks=1            # Run a single task
#SBATCH --cpus-per-task=12    # Number of CPU cores per task
#SBATCH --array=1-18          # Array range for 18 barcode folders
#SBATCH -p ei-medium           # Queue
#SBATCH --mem=80G            # Memory per CPU
#SBATCH --time=1-01:00:00       # Adjust runtime as needed
#SBATCH -o /ei/projects/1/17339e77-1a5d-440a-abde-70e615a0f98b/data/logs/process_fq_%A_%a.out # STDOUT
#SBATCH -e /ei/projects/1/17339e77-1a5d-440a-abde-70e615a0f98b/data/logs/process_fq_%A_%a.err # STDERR

cpus=$SLURM_CPUS_PER_TASK
source package 19a709c9-3ccf-4c82-acb9-0c0efeefff69  # restrander
source package 1041444f-cd25-4107-a5c7-5e86cb1728fe  # seqkit 2.8.2
restrander_config="/ei/projects/1/17339e77-1a5d-440a-abde-70e615a0f98b/data/results/20241031_IPSCs/references/PCB114.json"

# Define paths
# Accept configurable paths via command-line arguments or environment variables
base_dir="/ei/projects/1/17339e77-1a5d-440a-abde-70e615a0f98b/data/raw/20241031_IPSCs/_1/20241031_1631_P2S-00660-A_PAY86305_d9cfe10b/fastq_pass"
outdir="/ei/projects/1/17339e77-1a5d-440a-abde-70e615a0f98b/data/results/20241031_IPSCs/processed_fastq"
restrander_config="/ei/projects/1/17339e77-1a5d-440a-abde-70e615a0f98b/data/results/20241031_IPSCs/references/PCB114.json"

echo "Using base_dir: $base_dir"
echo "Using outdir: $outdir"
echo "Using restrander_config: $restrander_config"

# Create an array of barcode folder names
barcode_list=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18)
barcode=${barcode_list[$SLURM_ARRAY_TASK_ID-1]}
echo "Processing ${barcode}"
# Merge all fastq files in the barcode folder
mkdir -p $outdir/barcode$barcode
zcat "$base_dir/barcode$barcode"/*.fastq.gz | seqkit seq --remove-gaps --compress-level 5 --min-qual 8 -j $cpus \
-o "$outdir/barcode$barcode/merged.processed.fastq.gz"
seqkit stats -j $cpus --all -T "$outdir/barcode$barcode/merged.processed.fastq.gz" > \
"$outdir/barcode$barcode/merged.processed.stats.txt"

restrander "$outdir/barcode$barcode/merged.processed.fastq.gz" "$outdir/barcode$barcode/merged.processed.restranded.fastq.gz" "$restrander_config"

echo "Finished processing ${barcode}"