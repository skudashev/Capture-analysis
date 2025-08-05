# Merge all filtered BAMs into one
rule merge_filtered_bams:
    input:
        bams=expand(f"{OUTDIR}/alignments_2pass/{{sample}}.filt.bam", sample=SAMPLES)
    output:
        out_dir=directory(f"{OUTDIR}/tama"),
        merged=f"{OUTDIR}/tama/merged.secondpass.filt.bam"
    threads: config.get("threads", 4)
    shell: """
        mkdir -p $(dirname {output.merged})
        samtools merge -@ {threads} {output.merged} {input.bams}
        """

# Split merged BAM by genomic window
rule split_merged_bam:
    input:
        bam=f"{OUTDIR}/tama/merged.secondpass.filt.bam"
    output:
        touch(f"{OUTDIR}/tama/merged_bam_split.done")
    params:
        prefix=f"{OUTDIR}/tama/merged.secondpass.filt",
        script="scripts/split_bam_by_windows.py",
        window_size=config.get("tama_window_size", 10000000)
    shell: """
        python {params.script} {input.bam} {params.prefix} {params.window_size}
        touch {output}
        """

# Auto-detect region wildcards
TAMA_REGIONS = glob_wildcards(f"{OUTDIR}/tama/merged.secondpass.filt_{{region}}.sam").region

# TAMA collapse per split region
rule tama_collapse:
    input:
        sam=f"{OUTDIR}/tama/merged.secondpass.filt_{{region}}.sam"
    output:
        bed=f"{OUTDIR}/tama_collapse/merged.secondpass.filt_{{region}}_collapsed.bed"
    params:
        genome=config["genome"],
        prefix=f"{OUTDIR}/tama_collapse/merged.secondpass.filt_{{region}}"
    shell: """
        tama_collapse.py \
            -s {input.sam} \
            -f {params.genome} \
            -p {params.prefix} \
            -x no_cap -c 90 -i 80 \
            -icm ident_cov -log log_off -a 200 -z 100 -m 5 \
            -sjt 10 -lde 5 -e common_ends -rm low_mem
        """
