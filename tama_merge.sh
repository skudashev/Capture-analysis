#!/bin/bash
#SBATCH --job-name=tama_merge
#SBATCH --mail-type=ALL
#SBATCH --mail-user=youremail@example.com
#SBATCH -p ei-short
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=50g
#SBATCH --time=12:00:00
#SBATCH --output=logs/tama_merge_%A.log

source package 3f3c222b-852c-43b8-bcb4-8822c5d67fc2 # tama-1.0
tama_merge.py -f formatted_filelist.txt -p tama_collapse -e common_ends -a 20 -m 3 -z 20 -d merge_dup 

probes_gtf="/path/to/all_probes.gtf"
source package 4028d6e4-21a8-45ec-8545-90e4ed7e1a64 # bedtools -2.30.0
bedtools intersect -u -a tama_collapse.bed -b $probes_gtf > tama_transcriptome_probes_filt.bed

python tama_convert_bed_gtf_ensembl_no_cds.py tama_transcriptome_probes_filt.bed tama_iq_transcriptome.gtf