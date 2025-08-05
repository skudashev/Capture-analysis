rule tama_merge_and_convert:
    input:
        gtf=f"{OUTDIR}/tama/merged_capture.transcript_models.gtf",
        probes_gtf=config["capture_probes_gtf"],
        filelist=f"{OUTDIR}/tama/formatted_filelist.txt"
    output:
        filtered_gtf=f"{OUTDIR}/tama/merged_capture.transcript_models.filtered.gtf",
        bed=f"{OUTDIR}/tama/isoquant.bed",
        merged=f"{OUTDIR}/tama/tama_collapse.bed",
        probes_filt_bed=f"{OUTDIR}/tama/tama_transcriptome_probes_filt.bed",
        final_gtf=f"{OUTDIR}/tama/tama_iq_transcriptome.gtf"
    params:
        conda_env=config["tama_env"],
        tama_tools="scripts/tama_tools",  
    shell: """
        conda activate {params.conda_env}

        # Filter stranded entries
        awk '$7 == "+" || $7 == "-" {{print}}' {input.gtf} > {output.filtered_gtf}

        # Convert GTF to BED12
        python {params.tama_tools}/tama_format_gtf_to_bed12_ensembl.py \
            {output.filtered_gtf} {output.bed}

        # Run TAMA merge
        tama_merge.py -f {input.filelist} -p tama_collapse -e common_ends -a 20 -m 3 -z 20 -d merge_dup

        # Intersect with probes
        bedtools intersect -u -a {output.merged} -b {input.probes_gtf} > {output.probes_filt_bed}

        # Convert BED back to GTF for SQANTI
        python {params.tama_tools}/tama_convert_bed_gtf_ensembl_no_cds.py \
            {output.probes_filt_bed} {output.final_gtf}

        conda deactivate
    """
