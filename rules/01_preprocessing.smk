import os
import glob
import pandas as pd

configfile: 'config.yml'
samples = pd.read_csv(config["sample_sheet"], sep='\t', dtype=str).set_index(["sample_id"], drop=False)
SAMPLES = samples.index.tolist()
OUTDIR = config["outdir"]

rule merge_filter_fastq:
    input:
        lambda wildcards: sorted(glob.glob(os.path.join(samples.loc[wildcards.sample, "Path"], "*.fastq.gz")))
    output:
        fastq=f"{OUTDIR}/processed/{{sample}}/merged.processed.fastq.gz",
        stats=f"{OUTDIR}/processed/{{sample}}/merged.processed.stats.txt"
    threads: config.get("threads", 4)
    shell: """
        mkdir -p $(dirname {output.fastq})
        zcat {input} | seqkit seq --remove-gaps --compress-level 5 --min-qual 8 -j {threads} \
        -o {output.fastq}
        seqkit stats -j {threads} --all -T {output.fastq} > {output.stats}
        """

rule restrander:
    input:
        fastq=f"{OUTDIR}/processed/{{sample}}/merged.processed.fastq.gz"
    output:
        fastq=f"{OUTDIR}/processed/{{sample}}/merged.processed.restranded.fastq.gz"
    params:
        config=config["restrander_config"]
    shell: """
        mkdir -p $(dirname {output.fastq})
        restrander {input.fastq} {output.fastq} {params.config}
        """

rule map_to_genome:
    input:
        fastq=f"{OUTDIR}/processed/{{sample}}/merged.processed.restranded.fastq.gz"
    output:
        bam=f"{OUTDIR}/alignments/{{sample}}.sorted.bam",
        bai=f"{OUTDIR}/alignments/{{sample}}.sorted.bam.bai"
    params:
        genome=config["genome"]
    threads: config.get("threads", 4)
    shell: """
        mkdir -p $(dirname {output.bam})
        minimap2 -ax splice --cs=long -I 16G -G 500k -ub -k 14 --secondary=yes -t {threads} \
        {params.genome} {input.fastq} | samtools view -b | \
        samtools sort -O BAM -m 3G -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule score_junctions:
    input:
        bam=f"{OUTDIR}/alignments/{{sample}}.sorted.bam"
    output:
        bed=f"{OUTDIR}/junctions/{{sample}}.juncs.bed"
    params:
        genome=config["genome"],
        ref_juncs=config["ref_juncs"],
        sif=config["twopasstools"]
    threads: config.get("threads", 4)
    shell: """
        mkdir -p $(dirname {output.bed})
        singularity exec --env LC_ALL=C.UTF-8,LANG=C.UTF-8 {params.sif} \
        2passtools score -f {params.genome} -p {threads} -a {params.ref_juncs} \
        -m "GTAG|GCAG|ATAG|ATAC" -j 3 --keep-all-annot -o {output.bed} {input.bam}
        """

rule merge_junctions:
    input:
        beds=expand(f"{OUTDIR}/junctions/{{sample}}.juncs.bed", sample=SAMPLES)
    output:
        merged_bed=f"{OUTDIR}/junctions/merged_juncs.bed",
        filt_bed=f"{OUTDIR}/junctions/merged_juncs.filt.bed"
    params:
        genome=config["genome"],
        ref_juncs=config["ref_juncs"],
        twopasstools=config["twopasstools"],
        seed=42
    threads: config.get("threads", 4)
    log:
        f"{OUTDIR}/logs/merge_junctions.log"
    shell: """
        mkdir -p $(dirname {output.merged_bed})
        singularity exec --env LC_ALL=C.UTF-8,LANG=C.UTF-8 {params.twopasstools} \
        2passtools merge --keep-all-annot -f {params.genome} -p {threads} -a {params.ref_juncs} \
        -m "GTAG|GCAG|ATAG|ATAC" --classifier-type random_forest \
        --random-seed {params.seed} --verbosity DEBUG -o {output.merged_bed} {input.beds} >& {log}

        singularity exec --env LC_ALL=C.UTF-8,LANG=C.UTF-8 {params.twopasstools} \
        2passtools filter --exprs 'count > 2 and jad > 3 and (decision_tree_1_pred or decision_tree_2_pred)' \
        -o {output.filt_bed} {output.merged_bed}
        """

rule twopass_map_to_genome:
    input:
        fastq=f"{OUTDIR}/processed/{{sample}}/merged.processed.restranded.fastq.gz",
        filt_juncs=f"{OUTDIR}/junctions/merged_juncs.filt.bed"
    output:
        bam=f"{OUTDIR}/alignments_2pass/{{sample}}.sorted.bam",
        bai=f"{OUTDIR}/alignments_2pass/{{sample}}.sorted.bam.bai"
    params:
        genome=config["genome"],
        max_intron_size=config["max_intron_size"]
    threads: config.get("threads", 4)
    shell: """
        mkdir -p $(dirname {output.bam})
        minimap2 -ax splice -k 14 -t {threads} -ub -G {params.max_intron_size} \
        --junc-bed {input.filt_juncs} --junc-bonus 15 --secondary=yes \
        {params.genome} {input.fastq} | samtools sort -O BAM -m 3G -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule filter_read_ids:
    input:
        bam=f"{OUTDIR}/alignments_2pass/{{sample}}.sorted.bam"
    output:
        txt=f"{OUTDIR}/alignments_2pass/{{sample}}_reads_to_filter.txt"
    params:
        script="scripts/filter_RT.py",
        conda_env=config["lr_analysis_env"]
    shell: """
        conda activate {params.conda_env}
        python {params.script} --in_bam {input.bam} --out_txt {output.txt}
        conda deactivate
        """

rule filter_bam_reads:
    input:
        bam=f"{OUTDIR}/alignments_2pass/{{sample}}.sorted.bam",
        txt=f"{OUTDIR}/alignments_2pass/{{sample}}_reads_to_filter.txt"
    output:
        bam=f"{OUTDIR}/alignments_2pass/{{sample}}.filt.bam",
        bai=f"{OUTDIR}/alignments_2pass/{{sample}}.filt.bam.bai"
    threads: config.get("threads", 4)
    singularity:
        config["samtools"]
    shell: """
        mkdir -p $(dirname {output.bam})
        if [ ! -s {input.txt} ]; then
            echo "No reads to filter for sample {wildcards.sample}. Exiting."
            touch {output.bam}
            touch {output.bai}
        else
            samtools view -@ {threads} -bh -N ^{input.txt} -F 0x804 {input.bam} | \
            samtools sort -@ {threads} -o {output.bam}
            samtools index -@ {threads} {output.bam}
        fi
        """
