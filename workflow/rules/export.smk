rule export_amplicon_fasta:
    """One FASTA per assay of every valid amplicon (Detected and Primer Only),
    with flanking sequence, plus amplicon_index.csv. See export_amplicons.py."""
    input:
        records=expand(
            config["results_dir"] + "/amplicon_records/{accession}.tsv",
            accession=ACCESSIONS
        ),
        metadata=config["metadata"],
        assay_table=config["assay_table"],
    output:
        directory(config["results_dir"] + "/amplicon_fasta"),
    params:
        records_dir=config["results_dir"] + "/amplicon_records",
        group_by=" ".join(GROUP_BY),
    resources:
        mem_mb=4000,
    shell:
        """
        python {SCRIPTS}/export_amplicons.py \
            --records-dir "{params.records_dir}" \
            --metadata "{input.metadata}" \
            --assay-table "{input.assay_table}" \
            --group-by {params.group_by} \
            --out-dir "{output}"
        """
