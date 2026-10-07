# Per-amplicon records feed the amplicon FASTA export (rules/export.smk). They
# are intermediate (temp) and only produced when amplicon_fasta is on.
DETECT_OUTPUTS = {"detection": config["results_dir"] + "/amplicons/{accession}.csv"}
if AMPLICON_FASTA:
    DETECT_OUTPUTS["records"] = temp(
        config["results_dir"] + "/amplicon_records/{accession}.tsv")


rule run_detection:
    input:
        # BLAST_DIR is this run's own blast/ directory, or a previous run's
        # retained blast/ when re-scoring (see the Snakefile header).
        blast=BLAST_DIR + "/{accession}.tsv.gz",
        # Still required when re-scoring: detect.py re-reads the genome to
        # reconstruct BLAST end-trimmed hits to full oligo length.
        fna=config["assembly_dir"] + "/{accession}.fna",
        assay_table=config["assay_table"],
    output:
        **DETECT_OUTPUTS,
    params:
        max_primer_mismatches=config["max_primer_mismatches"],
        prime3_exact_nt=config["prime3_exact_nt"],
        max_probe_mismatches=config["max_probe_mismatches"],
        max_amplicon_size=config["max_amplicon_size"],
        store_amplicon_sequences=config["store_amplicon_sequences"],
        flank_bp=AMPLICON_FLANK_BP,
        records_arg=lambda wildcards, output: (
            f'--records-out "{output.records}"' if AMPLICON_FASTA else ""),
    resources:
        mem_mb=2000,
    shell:
        """
        python {SCRIPTS}/detect.py \
            --blast "{input.blast}" \
            --fna "{input.fna}" \
            --assay-table "{input.assay_table}" \
            --max-primer-mismatches {params.max_primer_mismatches} \
            --prime3-exact-nt {params.prime3_exact_nt} \
            --max-probe-mismatches {params.max_probe_mismatches} \
            --max-amplicon-size {params.max_amplicon_size} \
            --store-amplicon-sequences {params.store_amplicon_sequences} \
            --detection-out "{output.detection}" \
            --flank-bp {params.flank_bp} \
            {params.records_arg}
        """
