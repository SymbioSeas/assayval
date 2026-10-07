# Encoding-safe input tables + per-assay amplicon FASTA — design

Date: 2026-10-07. Status: implemented (user confirmed all four decision points). Header field `call=` shipped as `assembly_call=` for clarity alongside the per-copy `probe=`; CLI de-duplicates repeated Python warnings.

## Part 1 — encoding-safe user tables

**Problem.** `assay-val init` writes an ASCII `assay_table.csv`, but Excel re-saves
"CSV" in the platform codepage (Mac Roman on macOS, Windows-1252 on Windows).
Every reader used `encoding='utf-8-sig'`, so any non-ASCII free text (en dash,
accented author names, `°`) crashed the run with `UnicodeDecodeError`.

**Shared loader** `workflow/scripts/table_io.py`:

- `decode_table_bytes(data) -> (text, encoding)`: UTF-8 (BOM stripped) →
  UTF-16 (BOM only) → otherwise decode as both cp1252 and mac_roman and pick
  the higher plausibility score (non-ASCII char inside a word should be a
  letter; elsewhere it should be typographic punctuation/symbol). latin-1 as
  last resort. Non-UTF-8 → `warnings` line on stderr naming the guess and the
  "CSV UTF-8" fix.
- `read_table_text(path)` → decoded text; `sniff_delimiter(header)` accepts
  `,` `;` `\t`.
- `read_csv_rows(path)` → list of dicts; header names and cell values stripped
  (incl. NBSP U+00A0).
- `load_assay_table(path)` → validated rows: required columns
  `assay,fwd,rev,probe`; assay names non-empty, unique, no whitespace;
  `fwd`/`rev` non-empty; after modification-stripping every oligo is IUPAC
  (`ACGTRYSWKMBDHVN`, case-insensitive). Errors are `ValueError` naming row,
  assay and column.
- `read_metadata(path)` → `pandas.DataFrame` via the same decoder/delimiter.

Call sites switched: `prepare_oligos.load_assay_table`, `detect.run_detection`,
`summarize.load_assay_targets`, `summarize.load_assay_context`, summarize
metadata read. Fixes two latent silent failures: whitespace in an assay name
(BLAST truncates qseqid at the space → every hit unmatched) and trailing
whitespace in a primer (counted as a base → fails 3'-exact).

**Template side.** `assay-val init` writes `assay_table.csv` with a UTF-8 BOM
(Excel opens it as "CSV UTF-8" and a plain Save keeps that format).
`config/config.yaml` made pure ASCII. CLI reads config as UTF-8.

## Part 2 — per-assay amplicon FASTA

**Config** (defaults applied when absent):
`amplicon_fasta: true`, `amplicon_flank_bp: 50`.

**Flow.** `detect.py` (already holds the genome in memory) writes
`<run>/amplicon_records/{accession}.tsv` — one row per valid amplicon, `temp()`
— when `amplicon_fasta` is on. New rule `export_amplicon_fasta` →
`workflow/scripts/export_amplicons.py` reads all records + metadata + assay
table and writes `directory(<run>/amplicon_fasta)` containing
`{assay}_amplicons.fasta` for **every** assay (empty file when no amplicons) and
`amplicon_index.csv`.

**Which amplicons.** Every valid amplicon (fwd+rev pass mismatch / 3'-exact /
size filters): Detected and Primer Only alike.

**Sequence.** Oriented fwd-primer → rev-primer (minus-strand amplicons
reverse-complemented). Flanks `amplicon_flank_bp` each side, lowercase;
primer-to-primer amplicon uppercase; flanks clipped at contig ends (actual
lengths reported). Primer sites are the genome's sequence.

**Numbering.** `amp{k}of{n}`; k follows the order of `amplicon_starts` /
`contig_ids` in `amplicons/{accession}.csv`. Per-amplicon probe status.

**Header.**
`>{accession}_amp{k}of{n}_{Organism_sanitized} organism="{organism}" assay=…
call=… probe=yes|no|none contig={id}:{start}-{end}({strand}) amplicon_bp=…
flank_bp={up},{down} fwd_mm=… rev_mm=… probe_mm=…`
Sanitize: chars outside `[A-Za-z0-9.-]` → `_`, runs collapsed. Organism:
`organism_name` → first `group_by` column → `unknown`. File name: assay
sanitized the same way (collision → error).

**Also.** Manifest records both keys; CLI footer lists `amplicon_fasta/`;
README Outputs + config docs. `store_amplicon_sequences` unchanged.
