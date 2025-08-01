#!/bin/bash
#SBATCH --job-name=tama_collapse_array   
#SBATCH --mail-type=ALL       
#SBATCH --mail-user=youremail@example.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH -p ei-medium
#SBATCH --mem=40G
#SBATCH --array=1-10
#SBATCH --time=7-00:00:00
#SBATCH -o logs/tama_%A_%a.out
#SBATCH -e logs/tama_%A_%a.err

source ~/.bashrc
source package 3f3c222b-852c-43b8-bcb4-8822c5d67fc2  # tama-1.0

genome="/path/to/hg38_sequin.fa"
sam_dir="/path/to/split_sam_files"
out_dir="$sam_dir/tama_collapse"
mkdir -p $out_dir

sam_files=($sam_dir/*.sam)
sam_file=${sam_files[$SLURM_ARRAY_TASK_ID]}
chr=$(basename $sam_file .sam)

echo "Running TAMA collapse on ${sam_file}"
tama_collapse.py \
    -s ${sam_file} \
    -f $genome \
    -p ${out_dir}/${chr} \
    -x no_cap -c 90 -i 80 \
    -icm ident_cov -log log_off -a 200 -z 100 -m 5 \
    -sjt 10 -lde 5 \
    -e common_ends -rm low_mem
echo "Finished ${chr}"