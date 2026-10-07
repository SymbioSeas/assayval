# Encoding-safe tables + amplicon FASTA — implementation plan

Spec: `docs/superpowers/specs/2026-10-07-encoding-and-amplicon-fasta-design.md`.
TDD throughout: failing test → watch fail → minimal code → green.
Test runner: `PATH=~/miniforge3/envs/assayval/bin:$PATH python -m pytest tests/ -q`.

1. **table_io decode** — `tests/test_table_io.py`: UTF-8, UTF-8-BOM, UTF-16-BOM,
   mac_roman and cp1252 encodings of the real strings (`14–15`, `Fernández-No`,
   `60 °C`, `Caamaño-Antelo`) all decode to the original text; non-UTF-8 warns.
2. **table_io rows + delimiter** — `,` `;` `\t`; header/cell whitespace + NBSP
   stripped; CRLF fine.
3. **load_assay_table validation** — missing column, empty/duplicate/whitespace
   assay name, empty fwd, non-IUPAC oligo (after mod strip) each raise with the
   assay name; valid table with mods/degenerate bases passes.
4. **Switch call sites** — prepare_oligos, detect, summarize (assay table +
   metadata) use table_io; existing tests stay green; add a Mac-Roman
   end-to-end test through `prepare_oligos.load_assay_table` and
   `summarize.load_assay_context`.
5. **Templates** — init writes BOM; config.yaml pure ASCII (test both);
   CLI `_load_config` UTF-8.
6. **detect: strand + flanked extraction** — `find_valid_amplicons` records
   `strand`; `extract_flanked_amplicon(contigs, contig, start, end, strand,
   flank)` → (seq, up, down): plus/minus strand, contig-edge clipping, case.
7. **detect: amplicon records** — `build_amplicon_records(...)` rows incl.
   Primer Only, per-amplicon probe status, k/n order = detection CSV order;
   `--records-out/--flank-bp` CLI args.
8. **export_amplicons.py** — header format, sanitizing, organism fallback,
   empty FASTA for zero-amplicon assays, index CSV, filename collision error.
9. **Workflow wiring** — config keys (`_as_bool`, defaults), conditional
   detect output (temp), new rule, `rule all`, manifest params, CLI footer;
   Snakemake dry-run with flag on/off.
10. **Docs** — README (assay table/Excel note, config, Outputs), spec status.
11. **Verify** — full suite; real ~20-assembly Bacillales run with the user's
    assay table re-encoded to Mac Roman.
