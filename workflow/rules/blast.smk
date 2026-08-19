rule run_blast:
    input:
        oligos="resources/oligos/all_oligos.fasta",
        db_nhr="resources/blast_db/{accession}/{accession}.nhr",
    output:
        tsv=blast_output(config["results_dir"] + "/blast/{accession}.tsv.gz"),
    params:
        db="resources/blast_db/{accession}/{accession}",
        evalue=config["blast_evalue"],
        perc_identity=config["blast_perc_identity"],
        word_size=config["blast_word_size"],
        max_target_seqs=config.get("blast_max_target_seqs", 50000),
    resources:
        mem_mb=4000,
    log:
        config["results_dir"] + "/logs/{accession}.blast.log",
    shell:
        # Only the eight columns the detection engine actually consumes are
        # requested. In particular qseq/sseq are omitted: detect.py overwrites
        # the subject sequence anyway, re-extracting each hit's full-length
        # window from the assembly (see reconstruct_full_hits), so asking BLAST
        # for it doubles the size of every row to no purpose. Output is gzipped
        # because this is the dominant disk cost of a run and the format
        # compresses roughly ten-fold; pandas decompresses transparently.
        """
        set -o pipefail
        blastn \
            -task blastn-short \
            -query "{input.oligos}" \
            -db "{params.db}" \
            -evalue {params.evalue} \
            -perc_identity {params.perc_identity} \
            -word_size {params.word_size} \
            -max_target_seqs {params.max_target_seqs} \
            -strand both \
            -outfmt "6 qseqid sseqid length gapopen qstart qend sstart send" \
            2> "{log}" \
        | gzip -c > "{output.tsv}"
        """
