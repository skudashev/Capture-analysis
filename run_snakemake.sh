#!/bin/bash
#SBATCH --job-name=snake   # Job name
#SBATCH --mail-type=ALL             # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=sofia.kudasheva@earlham.ac.uk  # Where to send mail
#SBATCH -p ei-long # queue
#SBATCH -N 1 # number of nodes
#SBATCH -c 1 # number of cores
#SBATCH -t 10-00:00 # time (D-HH:MM)
#SBATCH --mem=10gb                     # Job memory request
#SBATCH -o /hpc-home/kudashev/logs/snakemake.%N.%j.out # STDOUT
#SBATCH -e /hpc-home/kudashev/logs/snakemake.%N.%j.err # STDERR

# source package 222eac79-310f-4d4b-8e1c-0cece4150333 #minimap2-2.24-(r1122)
source package /tgac/software/testing/bin/minimap2-2.24 
source package 638df626-d658-40aa-80e5-14a275b7464b #samtools-1.15.1
source package db01e69e-ca62-416d-a9bc-b28cb8139657 #sqanti3-5.1
source package 7775013e-6118-4b30-86dd-0d7056a466fd #PicardTools - 2.23.9
source package 19a709c9-3ccf-4c82-acb9-0c0efeefff69 #restrander v1.0.0
source package 4028d6e4-21a8-45ec-8545-90e4ed7e1a64 #bedtools-2.30.0
source snakemake-8.14.0_CBG #snakemake-8.14.0 - Gemy
source seqkit-2.5.1_CBG #seqkit-2.5.1 - Gemy

# source /usr/users/EI_ga011/kudashev/.bashrc
source /hpc-home/kudashev/.bashrc

cd /hpc-home/kudashev/pipeline
# snakemake  --profile profiles/smk_slurm_config --touch 
snakemake --profile  ./smk_slurm_config --keep-incomplete --rerun-triggers mtime