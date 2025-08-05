rule filter_gtf_stranded:
    input:
        gtf=f"{OUTDIR}/isoquant_output/OUT.transcript_models.gtf"
    output:
        gtf=f"{OUTDIR}/isoquant_output/OUT.transcript_models.stranded.gtf"
    shell: """
        if [ -f {output.gtf} ]; then
            echo "Filtered GTF file already exists: {output.gtf}"
        else
            mkdir -p $(dirname {output.gtf})
            awk '$7 == "+" || $7 == "-" {{print}}' {input.gtf} > {output.gtf}
            echo "Filtered GTF file created: {output.gtf}"
        fi
        """


rule sqanti3_qc:
    input:
        gtf=rules.filter_gtf_stranded.output.gtf,
        ref_gtf=config["gtf"],
        genome=config["genome"]
    output:
        report_dir=directory(f"{OUTDIR}/sqanti3_output"),
        isoforms=f"{OUTDIR}/sqanti3_output/OUT.transcript_models.stranded_corrected.fasta",
        gtf=f"{OUTDIR}/sqanti3_output/OUT.transcript_models.stranded_corrected.gtf",
        classification=f"{OUTDIR}/sqanti3_output/OUT.transcript_models.stranded_classification.txt",
        report=f"{OUTDIR}/sqanti3_output/OUT.transcript_models.stranded_SQANTI3_report.pdf"
    params:
        sr_bam="short_reads/bam_list.txt",
        sr_juncs="short_reads/junctions/",
        polyA=config["polyA"],
        polyA_motif=config["polyA_motif"],
        cage=config["CAGE"]
    threads: config.get("threads", 8)
    shell: """
        mkdir -p {output.report_dir}
        sqanti3_qc.py --cpus {threads} --report pdf --force_id_ignore --aligner_choice minimap2 \
            --polyA_motif_list {params.polyA_motif} \
            --CAGE_peak {params.cage} --isoAnnotLite \
            --min_ref_len 200 \
            --SR_bam {params.sr_bam} \
            --ratio_TSS_metric max \
            --polyA_peak {params.polyA} \
            --sites GTAG,GCAG,ATAC --dir {output.report_dir} \
            --coverage {params.sr_juncs} \
            {input.gtf} {input.ref_gtf} {input.genome}
        """


rule sqanti3_filter:
    input:
        isoforms=rules.sqanti3_qc.output.isoforms,
        gtf=rules.sqanti3_qc.output.gtf,
        classification=rules.sqanti3_qc.output.classification
    output:
        filtered_dir=directory(f"{OUTDIR}/sqanti3_filter")
    params:
        rules_json=config["sqanti_rules"]
    shell: """
        mkdir -p {output.filtered_dir}
        sqanti3_filter.py rules -j {params.rules_json} \
            --isoforms {input.isoforms} \
            --gtf {input.gtf} \
            --dir {output.filtered_dir} \
            {input.classification}
        """
