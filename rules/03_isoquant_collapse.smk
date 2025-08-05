# Generate labels.txt from the sample sheet for IsoQuant
rule isoquant_labels:
    input:
        sample_sheet=config["sample_sheet"]
    output:
        labels=f"{OUTDIR}/isoquant/labels.txt"
    params:
        script="scripts/make_labels.py",
        input_dir="alignments_2pass"
    shell:"""
        python {params.script} --sample_sheet {input.sample_sheet} \
               --input_dir {params.input_dir} --output_labels {output.labels}
        """

# Run IsoQuant once across all samples using the filtered BAMs and generated labels
rule isoquant_collapse:
    input:
        bams=expand(f"{OUTDIR}/alignments_2pass/{{sample}}.filt.bam", sample=SAMPLES),
        labels=f"{OUTDIR}/isoquant/labels.txt",
    output:
        log=f"{OUTDIR}/isoquant/isoquant.log",
    params:
        genome=config["genome"],
        genedb=config["genedb"],
        juncs=config["ref_juncs"],
        output_dir="isoquant_output"
    threads: config.get("threads", 8)
    shell:"""
        if [[ -f {output.log} ]]; then
            last_line=$(tail -n 1 {output.log})
            if [[ "$last_line" == *"IsoQuant pipeline finished"* ]]; then
                echo "IsoQuant pipeline already finished. Exiting."
                exit 0
            else
                echo "Resuming IsoQuant pipeline..."
                isoquant.py -t {threads} --output {params.output_dir} --resume
            fi
        else
            echo "Starting IsoQuant pipeline from scratch..."
            isoquant.py -t {threads} --data_type nanopore \
                --reference {params.genome} \
                --genedb {params.genedb} --complete_genedb \
                --bam {input.bams} --labels {input.labels} \
                --junc_bed_file {params.juncs} \
                --splice_correction_strategy conservative_ont \
                --model_construction_strategy sensitive_ont \
                --matching_strategy default \
                --polya_requirement auto --report_canonical all \
                --report_novel_unspliced false \
                --output {params.output_dir}
        fi
        """
