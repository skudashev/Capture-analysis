# authorship
__author__ = "Sofia Kudasheva"
__email__ = "sofia.kudasheva@earlham.ac.uk"

# Import necessary libraries
import os
from os import path
import pandas as pd

# Configurations
configfile: 'config.yml'
WORKDIR = os.getcwd()  # path to working directory
OUTDIR = config["outdir"]  # path to output directory
if not os.path.exists(OUTDIR):
    os.makedirs(OUTDIR)
logs_dir = path.join(WORKDIR,"logs")  
samples = pd.read_csv(config["sample_sheet"], sep='\t', dtype=str).set_index(["sample_id"], drop=False)
SAMPLES = samples.index.tolist()

# Rule includes
include: "rules/01_preprocessing.smk"
include: "rules/02_tama_collapse.smk"
include: "rules/03_isoquant_collapse.smk"
include: "rules/04_sqanti3.smk"

# Rule all to create output of the pipeline
rule all:
    input:
        # 01_preprocessing: filtered & restranded fastqs
        expand(f"{OUTDIR}/processed/{{sample}}/merged.processed.fastq.gz", sample=SAMPLES),
        # 01_preprocessing: filtered BAMs
        expand(f"{OUTDIR}/alignments_2pass/{{sample}}.filt.bam", sample=SAMPLES),
        # 02_tama_collapse: TAMA collapsed BEDs
        expand(f"{OUTDIR}/tama_collapse/merged.secondpass.filt_{{region}}_collapsed.bed", region=TAMA_REGIONS),
        f"{OUTDIR}/isoquant/isoquant.log"
        
        