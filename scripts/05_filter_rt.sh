#!/bin/bash
#SBATCH --job-name=filterRT   # Job name
#SBATCH --mail-type=ALL       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sofia.kudasheva@earlham.ac.uk  # Where to send mail
#SBATCH --ntasks=1            # Run a single task
#SBATCH --cpus-per-task=10    # Number of CPU cores per task
#SBATCH --array=1,5-18         # Array range for 18 barcode folders
#SBATCH -p ei-short           # Queue
#SBATCH --mem=25G            # Memory per CPU
#SBATCH --time=00:45:00       # Adjust runtime as needed
#SBATCH -o /ei/projects/1/17339e77-1a5d-440a-abde-70e615a0f98b/data/logs/filterRT_%A_%a.out # STDOUT
#SBATCH -e /ei/projects/1/17339e77-1a5d-440a-abde-70e615a0f98b/data/logs/filterRT_%A_%a.err # STDERR

# Define variables
scripts="/ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/scripts"
samtools="/ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/tools/samtools/samtools-1.20.img"
# Define the base directory for the barcode folders
output_dir="/ei/projects/1/17339e77-1a5d-440a-abde-70e615a0f98b/data/results/20241031_IPSCs/alignments/merged_barcodes"

# Create an array of barcode folder names
barcode_list=(01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18)
barcode=${barcode_list[$SLURM_ARRAY_TASK_ID-1]}
in_bam="${output_dir}/barcode${barcode}.sorted.bam"
out_txt="${output_dir}/barcode${barcode}_reads_to_filter.txt"

cpus=$SLURM_CPUS_PER_TASK

source ~/.bashrc
conda activate lr_analysis
python $scripts/filter_RT.py --in_bam $in_bam --out_txt $out_txt
conda deactivate

### Check if the reads_to_filter.txt file exists and is not empty

if [ ! -s $out_txt ]; then
    echo "No reads to filter for barcode $barcode. Exiting."
    exit 1
else
    singularity exec $samtools samtools view -@ $cpus -bh -N ^$out_txt -F 0x804 $in_bam | \
    singularity exec $samtools samtools sort -@ $cpus -o ${output_dir}/barcode${barcode}.filt.bam
    singularity exec $samtools samtools index -@ $cpus ${output_dir}/barcode${barcode}.filt.bam
fi

