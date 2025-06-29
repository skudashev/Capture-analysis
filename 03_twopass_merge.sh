#!/bin/bash
#SBATCH --job-name=twopass   # Job name
#SBATCH --mail-type=ALL       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sofia.kudasheva@earlham.ac.uk  # Where to send mail
#SBATCH --ntasks=1            # Run a single task
#SBATCH --cpus-per-task=12    # Number of CPU cores per task
#SBATCH -p ei-medium           # Queue
#SBATCH --mem=200G            # Memory per CPU
#SBATCH --time=1-00:45:00       # Adjust runtime as needed
#SBATCH -o /ei/projects/1/17339e77-1a5d-440a-abde-70e615a0f98b/data/logs/twopass_%A.out # STDOUT
#SBATCH -e /ei/projects/1/17339e77-1a5d-440a-abde-70e615a0f98b/data/logs/twopass_%A.err # STDERR

source package 222eac79-310f-4d4b-8e1c-0cece4150333  #minimap2 2.24
source package 638df626-d658-40aa-80e5-14a275b7464b #samtools-1.15.1

# Define the base directory for the barcode folders
output_dir="/ei/projects/1/17339e77-1a5d-440a-abde-70e615a0f98b/data/results/20241031_IPSCs/alignments"
genome="/ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/references/refdata-gex-GRCh38-2024-A/fasta/genome.fa"
ref_juncs="/ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/references/gencode.v46.annotation.annotated_juncs.bed"
twopasstools="/ei/projects/3/3e4d9f24-546b-4e65-b3bd-e14aa572cf87/data/SK/annotation/tools/2passtools/2passtools.sif"
cpus=$SLURM_CPUS_PER_TASK
random_seed=42
juncs_list=$(ls ${output_dir}/*/*barcode*juncs.bed)

singularity exec --env LC_ALL=C.UTF-8,LANG=C.UTF-8 $twopasstools 2passtools merge --keep-all-annot -f $genome -p $cpus -a $ref_juncs -m "GTAG|GCAG|ATAG|ATAC" \
        --classifier-type random_forest --random-seed $random_seed --verbosity DEBUG -o ${output_dir}/merged_juncs.bed $juncs_list

singularity exec --env LC_ALL=C.UTF-8,LANG=C.UTF-8 $twopasstools 2passtools filter --exprs 'count > 2 and jad > 3 and (decision_tree_1_pred or decision_tree_2_pred)' \
        -o ${output_dir}/merged_juncs.filt.bed ${output_dir}/merged_juncs.bed