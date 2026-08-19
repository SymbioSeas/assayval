# assay-design: clade-specific ortholog discovery

Companion tool to [AssayVal](../README.md). `assay-design` parses a
[Panaroo](https://gtonkinhill.github.io/panaroo/) pangenome to find orthologs
that are **conserved within** a clade and **specific to** the clade (absent from the
rest of the dataset), then extracts representative protein and nucleotide
sequences. These clade-specific genes are the candidate targets from which the
Vpop dPCR assays were designed; the resulting amplicons are then validated with
AssayVal.

## Workflow

```
Panaroo output ──▶ 
Stage 1: identify conserved & specific orthologs (parse_group_specific_orthologs.py) ──▶ 
Stage 2: extract representative sequences (extract_proteins.py)
```

Both stages are wrapped by the `assay-design` command.

## Environment

Uses AssayVal's single top-level environment, so no separate install is required:

```bash
cd assayval
conda env create -f environment.yaml   # if not already created
conda activate assayval
pip install -e .                       # if not already installed
```

Installing the package puts both the `assay-val` and `assay-design` commands on
your PATH. Only `python`, `pandas`, and `numpy` are needed for this tool (all
included in that environment).

## Inputs

| Input                                     | Format                                                                                                                                                                                                                                                       |
| ----------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| **Presence/absence matrix** (`--matrix`)  | Panaroo `gene_presence_absence_roary` output (`.csv` or `.tsv`). Metadata columns (`Gene`, `Annotation`, …) are detected automatically; all other columns are treated as genomes. Each cell holds the locus ID for that genome (empty = gene absent).        |
| **Inclusion groups** (`--isolates-dir`)   | A directory of `.txt` files, one per group. Each file lists the genome names (one per line) that make up that inclusion group. Names must match matrix column headers. The filename stem names the group in all outputs.                                     |
| **Gene data** (`--gene-data`)             | Panaroo `gene_data.csv` with columns `gff_file`, `annotation_id`, `prot_sequence`, `dna_sequence`, `gene_name`, `description`. **All extracted sequences come from this file** — assay-design does not read separate per-genome FASTA (`.faa`) files. Only needed for Stage 2. |
| **Representatives** (`--representatives`) | Tab-separated file with a header row `group_stem` and `representative_assembly`, then one row per group. `group_stem` must match an inclusion-group filename stem. `representative_assembly` is a **genome name** — it must match a matrix column header **and** a `gene_data` `gff_file` value exactly (it is *not* a file path or `.faa` filename). Only needed for Stage 2. |

For a given group, the **exclusion set** is every genome column in the matrix
that is not in that group.

### How the identifiers connect

The same genome names appear in four places and must match exactly:

- column headers in the presence/absence matrix,
- `gff_file` values in `gene_data.csv`,
- lines in each inclusion-group `.txt` file, and
- the `representative_assembly` column of `representatives.tsv`.

So a `representatives.tsv` for a group `Lbrevis_Vjas` whose representative genome
is `LB14LO7` is:

```
group_stem	representative_assembly
Lbrevis_Vjas	LB14LO7
```

(tab-separated). Use the genome name `LB14LO7`, **not** a path like
`.../Lbrevis_faa/LB14LO7.faa`. Stage 2 then pulls that genome's protein and
nucleotide sequences for each ortholog directly from `gene_data.csv`.

## Usage

Run both stages on one dataset with the `assay-design` command (installed on your
PATH by the AssayVal conda environment):

```bash
assay-design \
    --matrix          gene_presence_absence.csv \
    --isolates-dir    isolate_groups \
    --gene-data       gene_data.csv \
    --representatives representatives.tsv \
    --output-dir      results
```

Useful flags:

- `--threshold 1.0` — strict all-or-nothing (present in *all* inclusion, absent from *all* exclusion)
- `--skip-stage2` — ortholog identification only (no `--gene-data`/`--representatives` needed)
- `--min-match-rate 0.5` — Stage 2 guard (see below)

Or run the stages individually:

```bash
python parse_group_specific_orthologs.py --matrix gene_presence_absence.csv --isolates-dir isolate_groups
python extract_proteins.py --orthologs-dir <dataset>_orthologs --gene-data gene_data.csv --representatives representatives.tsv
```

## Worked example

A tiny, ready-to-run example (subsampled from the FMS2026 dataset) lives in
[`example/`](example/):

```bash
cd example
bash run_example.sh
```

It runs in seconds and produces conserved/specific ortholog tables and
representative FASTAs for three *V. mediterranei* clades. See
[`example/README.md`](example/README.md).

## Outputs

Per inclusion group:

| File | Contents |
|------|----------|
| `<group>_conserved_orthologs.tsv` | Orthogroups conserved in the group, with inclusion/exclusion coverage metrics (`n_/frac_inclusion_present`, `n_/frac_exclusion_present`, `exclusion_assemblies_present`). |
| `<group>_specific_orthologs.tsv` | The conserved subset that is also specific to the group (clade-specific candidates). |
| `<group>_conserved_proteins.faa` | Representative protein sequences for the conserved orthologs (Stage 2). |
| `<group>_conserved_genes.fna` | Representative nucleotide sequences for the conserved orthologs (Stage 2). |
| `<group>_specific_proteins.faa` | Representative protein sequences for the group-specific orthologs — the clade-specific candidate targets (Stage 2). |
| `<group>_specific_genes.fna` | Representative nucleotide sequences for the group-specific orthologs (Stage 2). |

The `specific` FASTAs are the subset of the `conserved` FASTAs whose orthologs
are also absent from the exclusion genomes; they are the most direct
candidate-target set for assay design.

## Thresholds and rationale

**`--threshold` (default 0.9)** governs both calls:

- **Conserved** — present in **≥ 90%** of the group's inclusion genomes. Requiring
  near-universal (rather than strictly universal) presence tolerates draft-genome
  incompleteness and occasional annotation/clustering dropout. At the same time, it ensures the gene is a stable feature of the clade, a prerequisite for a sensitive assay.
- **Specific** — conserved **and** present in **≤ 10%** of the exclusion genomes
  (i.e. `1 − threshold`). Allowing a small exclusion fraction absorbs rare
  horizontal gene transfer and annotation noise without letting commonly shared
  genes through, a prerequisite for a specific assay.

Set `--threshold 1.0` for the strict definition (present in every inclusion
genome, absent from every exclusion genome). Values below ~0.8 are not
recommended: we find they admit orthologs too patchy to make reliable assay targets.

**`--min-match-rate` (default 0.5)** is a Stage 2 safety guard. Sequences are
extracted by matching each conserved ortholog's representative locus ID against
`gene_data` `annotation_id`s. If a genome's annotation-ID format in `gene_data`
differs from the presence/absence matrix, fewer than this fraction of IDs will
match and the group is **skipped with a warning** rather than producing a
misleadingly truncated FASTA. Raise it to be stricter about extraction
completeness.

## Note on functional annotation

Downstream functional annotation of the extracted proteins (e.g. with
[EggNOG-mapper](http://eggnog-mapper.embl.de/)) is outside the scope of this tool, but was useful for determining which
candidate gene targets to consider in developing the Vpop assay suite. The `.faa` files 
produced by Stage 2 are directly usable as input to EggNOG-mapper if desired.