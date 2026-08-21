rule aggregate_report:
    input:
        # ACCESSIONS is resolved once in the Snakefile: from assembly_dir on a
        # normal run, from the cached BLAST output when re-scoring.
        detections=expand(
            config["results_dir"] + "/amplicons/{accession}.csv",
            accession=ACCESSIONS
        ),
        metadata=config["metadata"],
        assay_table=config["assay_table"],
    output:
        matrices=expand(config["results_dir"] + "/reports/{gcol}_detection_matrix.csv", gcol=GROUP_COLS),
        heatmaps=expand(config["results_dir"] + "/reports/figures/{gcol}_detection_heatmap.pdf", gcol=GROUP_COLS),
        detection_long=config["results_dir"] + "/reports/detection_summary_long.csv",
        assay_summary_xlsx=config["results_dir"] + "/reports/assay_summary.xlsx",
        detection_by_assembly=config["results_dir"] + "/reports/detection_by_assembly.csv",
        assay_performance=config["results_dir"] + "/reports/assay_performance.csv",
        manifest=config["results_dir"] + "/reports/run_manifest.txt",
    params:
        amplicons_dir=config["results_dir"] + "/amplicons",
        reports_dir=config["results_dir"] + "/reports",
        group_by=" ".join(GROUP_BY),
        max_primer_mismatches=config["max_primer_mismatches"],
        prime3_exact_nt=config["prime3_exact_nt"],
        max_probe_mismatches=config["max_probe_mismatches"],
        max_amplicon_size=config["max_amplicon_size"],
        store_amplicon_sequences=config["store_amplicon_sequences"],
        keep_blast=KEEP_BLAST,
        keep_logs=KEEP_LOGS,
        # Recorded in run_manifest.txt so a re-scored run states which run's
        # BLAST output it was derived from. Empty on a normal run.
        rescored_from=RESCORE_FROM or "",
        # BLAST search provenance; empty on a re-scored run (see the Snakefile).
        blast_params=BLAST_PARAMS,
    resources:
        mem_mb=16000,
    shell:
        """
        python {SCRIPTS}/summarize.py \
            --amplicons-dir {params.amplicons_dir} \
            --metadata "{input.metadata}" \
            --assay-table "{input.assay_table}" \
            --reports-dir {params.reports_dir} \
            --group-by {params.group_by} \
            --max-primer-mismatches {params.max_primer_mismatches} \
            --prime3-exact-nt {params.prime3_exact_nt} \
            --max-probe-mismatches {params.max_probe_mismatches} \
            --max-amplicon-size {params.max_amplicon_size} \
            --store-amplicon-sequences {params.store_amplicon_sequences} \
            --keep-blast {params.keep_blast} \
            --keep-logs {params.keep_logs} \
            --rescored-from "{params.rescored_from}" \
            --blast-params "{params.blast_params}"
        """
