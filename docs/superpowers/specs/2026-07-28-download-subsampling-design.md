# Per-taxon random subsampling for `download-assemblies` — Design

**Status:** Approved (design), pending implementation plan.

## Goal

Let a user cap how many assemblies are downloaded per taxon, instead of always
downloading every assembly that matches the level/source filters. This supports
building manageable, representative assembly sets when downloading *all* genomes
for a taxon is unnecessary (e.g. specificity/background taxa for a multi-target
panel), while leaving the "download everything" behavior (e.g. the *Vibrionaceae*
use case) unchanged.

## Non-goals (deliberately out of scope)

- `--accessions <file>` curated-list input. Curated, allele-diverse positive
  sets (e.g. for AMR-gene sensitivity) are handled by side-loading `.fna` files
  into `assembly_dir` — AssayVal reads any `.fna`, so no code is needed. This is
  documented, not built.
- Quality- or representative-weighted selection (e.g. prefer complete/reference
  genomes). Selection is uniform-random within the already-applied
  `--assembly-level` / `--assembly-source` filters.
- Per-`-t` variable limits in a single command (`-t Klebsiella:20`). The compose
  pattern below covers this need without new CLI syntax.

## User-facing behavior

Two new flags on `download-assemblies` (`scripts/download/download_assemblies.sh`),
in the existing short-flag `getopts` style (`-t/-o/-l/-s/-e/-k`):

| Flag | Meaning | Default |
|------|---------|---------|
| `-n N` | Cap **each `-t` taxon** in this invocation at `N` randomly-selected assemblies. | unset → download all (unchanged) |
| `-r SEED` | Integer seed for the random selection. | `0` (reproducible by default) |

### Semantics

- **Per taxon, per invocation.** With `-t A -t B -n 20`, up to 20 are selected
  from A and 20 from B (each taxon sampled independently), *then* the existing
  cross-taxon de-duplication runs — so a genome shared between A and B is
  downloaded once (final count for an overlapping taxon may be `< N`).
- **Sampling is applied after** the `--assembly-level` and `--assembly-source`
  filters (it samples from the already-filtered candidate set).
- **`N ≥ available`** → take all available (no error).
- **`N ≤ 0` or non-integer** → validation error and exit non-zero, mirroring the
  existing `-s` value validation. `-r` requires an integer.
- **No `-n`** → every flag related to sampling is inert; behavior is byte-for-byte
  the current behavior.

### Compose pattern (covers the mixed-taxa use case)

Because the downloader is resume-aware and de-duplicates into one output
directory, different per-taxon N is achieved by composing invocations rather than
a more complex CLI:

```bash
download-assemblies -t Klebsiella          -n 20 -r 1 -o assemblies/
download-assemblies -t "Escherichia coli"  -n 20 -r 1 -o assemblies/
download-assemblies -t Vibrionaceae               -o assemblies/   # all
# -> assemblies/ holds 20 Klebsiella + 20 E. coli + all Vibrionaceae
```

## Implementation approach

### New module: `scripts/download/subsample.py`

A small, unit-testable Python helper (standard library only, matching
`parse_metadata.py`'s style). Core function:

```
sample_records(records: list[dict], n: int, seed: int) -> list[dict]
```

- Determinism: sort the input records by `accession` (canonical order), then
  apply a seeded `random.Random(seed).shuffle`, then take the first `n`.
  Sort-then-seeded-shuffle means a given input set + seed always yields the same
  subset, independent of the order NCBI returned them in.
- `n <= 0` or `n >= len(records)` → return all records (sorted, unshuffled order
  is fine; the contract is "the set," not the order).
- A thin CLI wrapper reads a JSON-Lines file (one `datasets summary` record per
  line), writes the sampled subset as JSON-Lines: e.g.
  `python3 subsample.py --jsonl <in> --n N --seed S --out <out>`.

### Wiring into `download_assemblies.sh`

Inside the existing per-taxon query loop
([download_assemblies.sh:139](../../../scripts/download/download_assemblies.sh)):

1. Write each taxon's `datasets summary` output to a per-taxon temp JSONL.
2. If `-n` is set, pass that temp file through `subsample.py` (with `-r`/seed)
   and append the sampled JSONL to `SUMMARY_ALL`; otherwise append it unchanged.
3. Log per taxon: `sampled N of M available (seed=S)`.

The existing global de-duplication + accession-list + download steps are
unchanged: they operate on `SUMMARY_ALL` exactly as today.

### Argument parsing / validation

- Extend `getopts` to `"t:o:l:s:n:r:e:k:h"`.
- Add `N`/`SEED` defaults (`N` empty = "all"; `SEED=0`).
- Validate: if `-n` is set, it must be a positive integer; if `-r` is set, it
  must be an integer. Emit a clear error + exit 1 on violation (same shape as the
  `-s` validation block).
- Document `-n`/`-r` in the header comment block that `usage()` prints, so they
  appear in `-h`.

## Reproducibility

The seed makes the selection deterministic **given the same NCBI result set**.
Because NCBI's holdings change over time, the authoritative record of what was
actually selected is the run's `metadata.csv` (one row per downloaded assembly)
plus `download.log` (which will now also record the seed and the sampled/available
counts per taxon). This matches how the pipeline already treats provenance.

## Testing

Unit tests for `sample_records` (no network):

- **Determinism:** same `(records, n, seed)` → identical selected accession set.
- **Seed sensitivity:** two different seeds over a sufficiently large input give
  (generally) different subsets.
- **`n >= len`:** returns all records.
- **`n <= 0`:** returns all records (or is rejected upstream — the CLI validates;
  the function's contract is "all" for non-positive to be safe).
- **Subset property:** every selected accession is present in the input.
- **Size:** exactly `min(n, len)` records returned for valid positive `n`.

CLI-level checks (mirroring the existing `-s` tests in `tests/test_download_cli.py`):

- `-h` lists `-n` and `-r`.
- A bad `-n` value (e.g. `-n abc` or `-n 0`) exits non-zero with a helpful message.

## Documentation

Add a short **"Scaling down large taxa"** subsection to the README's
[Downloading assemblies](../../../README.md) section:

- The compose pattern above (mixed per-taxon N via repeated invocations).
- The reproducibility note (seed + `metadata.csv`).
- A one-line pointer: for a *curated* set (e.g. allele-diverse AMR-gene carriers),
  download those accessions directly and drop the `.fna` files into
  `assembly_dir` — random subsampling is for scaling down broad taxa, not for
  curating a specific positive set.

## Files touched

| File | Change |
|------|--------|
| `scripts/download/subsample.py` | **New** — testable sampler + thin CLI. |
| `scripts/download/download_assemblies.sh` | Add `-n`/`-r`, validation, per-taxon sampling in the query loop, logging. |
| `tests/test_subsample.py` | **New** — unit tests for `sample_records`. |
| `tests/test_download_cli.py` | Add `-h` lists `-n`/`-r`; bad `-n` exits non-zero. |
| `README.md` | "Scaling down large taxa" subsection. |
