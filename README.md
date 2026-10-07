# AssayVal: in silico PCR assay validation

**AssayVal** evaluates the sensitivity and specificity of PCR assays (including probe-based assays designed for dPCR/qPCR) against a user-provided set of genome assemblies. For each assay, it identifies valid amplicons using BLAST-based primer alignment and reports detection calls, mismatch counts, and amplicon sizes per assembly.

## Two tools in this repository

| Tool                         | Location                                    | Purpose                                                                                                                                                                                                                                               |
| ---------------------------- | ------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **AssayVal**      | repository root (`workflow/`, `config/`, …) | Validate PCR assays in silico against a set of genome assemblies.                                                                                                                                                                                     |
| **assay-design** | [`assay-design/`](assay-design/)            | Identify clade-conserved / clade-specific orthologs from a Panaroo pangenome and extract representative sequences. This is the candidate-gene discovery step used to develop the Vpop assays. See [`assay-design/README.md`](assay-design/README.md). |

Both tools share the single conda environment defined in `environment.yaml`.

Once installed (below), each tool has a command on your PATH: **`assay-val`**
and **`assay-design`**. Run either with `--help`. Start a new AssayVal analysis
with `assay-val init` (see [Quick start](#quick-start)).

## Features

- Evaluates primer and probe binding across a user-provided set of genome assemblies
- Supports probe-based assays (hydrolysis probes, e.g., dPCR/qPCR) and probe-free assays (SYBR/dsDNA)
- IUPAC degenerate base support in all primer and probe sequences
- Handles primer binding on either strand of an assembly
- Configurable mismatch tolerances and 3′-exact match requirements
- Outputs per-assay detection summaries, amplicon details, and detection heatmaps
- Exports one FASTA per assay of every amplicon found (with flanking sequence), ready for probe/primer design in e.g. Geneious
- Accepts assay tables saved straight from Excel on Mac, Windows or Linux (legacy encodings and `;` delimiters)
- Runs locally or on SLURM HPC clusters via Snakemake profiles

## Requirements

- [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://mamba.readthedocs.io/) — all dependencies (including the `datasets` CLI) are installed by `environment.yaml`.
- **Platforms:** Linux and macOS natively; Windows via [WSL2](#installing-on-windows).

## System requirements

Disk and memory usage scale with the number and size of input assemblies. Estimates
below are from the manuscript's *Vibrionaceae* runs (RefSeq assemblies average
~5 MB each). Raw BLAST output is the dominant transient cost; with the default
`keep_blast: false` it is deleted as the run proceeds.

| Component | Per assembly | Draft set (10,715) | Complete set (927) |
|-----------|-------------|--------------------|--------------------|
| Input assemblies (`.fna`) | ~5 MB | 49 GB | 4.3 GB |
| Cached BLAST DBs (`resources/blast_db/`) | ~1.3 MB | 14 GB | ~1.2 GB |
| Raw BLAST output (transient) | ~15 MB | 158 GB | ~14 GB |
| Reports + amplicons | — | ~0.2 GB | small |
| **Peak disk (`keep_blast: false`)** | | **~65 GB** | **~6 GB** |
| **Peak disk (`keep_blast: true`)** | | **~220 GB** | **~20 GB** |
| **Peak RAM** | | **~16 GB** | **~16 GB** |

Rule of thumb: budget roughly **6–7 MB of transient disk per assembly** with
`keep_blast: false`, or ~21 MB per assembly if retaining raw BLAST, plus ~16 GB
RAM for the final aggregation step.

## Installation

```bash
git clone https://github.com/SymbioSeas/assayval.git
cd assayval
conda env create -f environment.yaml
conda activate assayval
pip install -e .
```

The pipeline runs inside this activated environment; the Snakemake profiles set
`use-conda: false` so no per-rule environments are built. The `pip install -e .`
step installs the package and puts three commands onto your PATH: **`assay-val`**, 
**`assay-design`**, and **`download-assemblies`** (the assembly-download helper).


### Installing on Windows

AssayVal runs on Windows through **WSL2** (Windows Subsystem for Linux), which
provides a real Linux environment. Native Windows is not supported because NCBI
BLAST+ is distributed for Linux/macOS only via Bioconda.

1. Install WSL2 and a Linux distribution (e.g. Ubuntu) — in an admin PowerShell:
   ```powershell
   wsl --install
   ```
   then reboot and finish the Ubuntu first-run setup. (See Microsoft's
   [WSL install guide](https://learn.microsoft.com/windows/wsl/install).)
2. **Inside the WSL/Ubuntu shell**, install Miniforge (conda):
   ```bash
   wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
   bash Miniforge3-Linux-x86_64.sh
   ```
3. From that same Linux shell, follow the standard [Installation](#installation)
   steps above. Everything else in this README then applies unchanged — you are
   running Linux.

Note: Keep your data and the repository **inside the WSL filesystem** (e.g. under
`~/`), not on the mounted Windows `/mnt/c/...` drive, for much faster BLAST I/O.


## Quick start

Each analysis lives in its own **analysis directory**: the folder that holds that
project's `config.yaml`, `assay_table.csv`, and `assemblies/`, and where results are
written. All commands below are run from inside it.

### 1. Create an analysis directory

Make a new directory for the project and drop the starter files into it with
`assay-val init`:

```bash
mkdir my_analysis && cd my_analysis
assay-val init
```

This writes two files you edit for your project in the following steps:

| File | Contents |
|------|----------|
| `config.yaml` | The default pipeline configuration (paths, detection thresholds, grouping, output options), identical to [`config/config.yaml`](config/config.yaml). |
| `assay_table.csv` | An assay table with all supported columns and two example rows (one probe-based, one probe-free) to replace with your assays. |

`assay-val init` never overwrites an existing `config.yaml` or `assay_table.csv`
(it skips them and exits non-zero); pass `--force` to replace them with fresh
copies. Use `--directory DIR` to set up a directory other than the current one
(it is created if missing).

### 2. Download assemblies

Use the `download-assemblies` command to download RefSeq assemblies for your taxon of interest:

```bash
download-assemblies -t "Vibrionaceae" -o assemblies/
```

This downloads all RefSeq assemblies (complete through contig level) for the specified taxon and writes a `metadata.csv` to the output directory. See [Downloading assemblies](#downloading-assemblies) for options and HPC usage.

### 3. Configure the pipeline

Edit the `config.yaml` created by `assay-val init`. Its defaults already match the
layout above (`assemblies/`, `assemblies/metadata.csv`, `assay_table.csv`), so
usually you only review the detection thresholds and grouping. All paths are
relative to the analysis directory:

```yaml
assembly_dir: "assemblies"          # directory containing .fna files
metadata: "assemblies/metadata.csv" # metadata CSV from download step
results_dir: "results"              # all outputs written here
assay_table: "assay_table.csv"      # your assay definitions

max_primer_mismatches: 2            # mismatches allowed per primer
prime3_exact_nt: 1                  # 3′-terminal bases that must match exactly
max_probe_mismatches: 1             # mismatches allowed in probe
max_amplicon_size: 500              # maximum amplicon size (bp)

# group_by: metadata column(s) used to group assemblies in the reports (a
# detection matrix + heatmap is written per column). Optional: leave unset to
# auto-group by NCBI ANI-derived species (requires metadata from
# download-assemblies). See "Metadata format" below.
# group_by: "species"
# group_by: ["species", "phenotype"]

amplicon_fasta: true                # one FASTA of all amplicons per assay
amplicon_flank_bp: 50               # flanking bp each side in those FASTAs
```

### 4. Prepare your assay table

Edit the `assay_table.csv` created by `assay-val init`: replace the two example rows
with your assays, one row per assay (see [Assay table format](#assay-table-format)).

### 5. Run

`assay-val` is installed on your PATH (see [Installation](#installation)). Run it
from the analysis directory, which now contains your `assemblies/` and the edited
`config.yaml` and `assay_table.csv`:

```bash
assay-val --run-name Vpop
```

Results are written to `<results_dir>/Vpop_<date>/` (`amplicons/`, `blast/`,
`reports/`, and `run.log`), where `<results_dir>` comes from `results_dir` in
`config.yaml` (default `results`). Options:

| Flag | Default | Description |
|------|---------|-------------|
| `--run-name NAME` | `results` | Names the results directory: `results/<NAME>_<date>/`. |
| `--directory DIR` | current dir | Analysis directory (assemblies, config, outputs). |
| `--configfile FILE` | `<dir>/config.yaml` | Pipeline configuration. |
| `--force` | off | Reuse today's `<NAME>_<date>/` and resume unfinished work. |
| `--cores N` | 8 | CPU cores. |
| `--set KEY=VALUE` | – | Override one config value for this run; repeatable. Mainly for detection thresholds, e.g. `--set max_primer_mismatches=1`. |
| `--rescore-from RUN_DIR` | – | Re-score a previous run's retained BLAST output instead of running BLAST again. See [Threshold sensitivity analysis](#threshold-sensitivity-analysis-re-scoring). |

Re-running the same name on the same day without `--force` creates
`<NAME>_<date>_2`, `_3`, … so previous results are never overwritten. Anything
after `--` is passed straight to Snakemake (e.g. `assay-val --run-name Vpop -- -n`
for a dry run). Set `keep_blast: true` in `config.yaml` to retain the raw
per-assembly BLAST output (see [System requirements](#system-requirements)), or
`keep_logs: true` to retain the per-assembly `logs/` directory (both are deleted
on success by default; `run.log` is always kept).

**Advanced (direct Snakemake / SLURM):** invoke the workflow directly, passing the
config explicitly:

```bash
snakemake --configfile config/config.yaml --profile workflow/profiles/local \
  --config results_dir=results/Vpop_manual
# SLURM: swap in --profile workflow/profiles/slurm
```

---

## Assay table format

The assay table is a CSV file with the following columns. `assay-val init` writes a
starter `assay_table.csv` with every column below already in place (the template is
[`config/assay_table_template.csv`](config/assay_table_template.csv)).

| Column  | Required | Description                                                                               |
| ------- | -------- | ----------------------------------------------------------------------------------------- |
| `assay` | Yes      | Unique assay name (used in all output files)                                              |
| `fwd`   | Yes      | Forward primer sequence (5′-3′)                                                           |
| `rev`   | Yes      | Reverse primer sequence (5′-3′, same orientation as fwd - AssayVal handles RC internally) |
| `probe` | Column yes, value no | Probe sequence (5′-3′). The **column must be present**, but leave the value empty to declare a probe-free (SYBR/dsDNA-dye) assay. See [Probe-free assays](#probe-free-assays). |
| `target_group` | No (column may be omitted) | The metadata group this assay is designed to detect, as `column:value` (e.g. `phenotype:protective`, `species:Vibrio mediterranei`). A bare value (no colon) is matched against the primary `group_by` column. Drives `assay_performance.csv`. |
| `target_gene` | No (column may be omitted) | Free-text gene/target label. Not used in detection; carried through to `assay_performance.csv` if present. |
| `max_primer_mismatches`, `prime3_exact_nt`, `max_probe_mismatches`, `max_amplicon_size` | No (columns may be omitted) | Per-assay overrides of the [detection thresholds](#detection-thresholds) in `config.yaml`. Leave a cell blank to use the config value; a whole number overrides it for that assay only. See [Per-assay thresholds](#per-assay-thresholds). |

**Sequence notation:** sequences are written 5′→3′ in any case, using IUPAC codes
(A C G T and the ambiguity codes R Y S W K M B D H V N), plus the modification
notation below. See [Oligo notation](#oligo-notation-lna-and-modifications) for the full
list.

**Which columns must exist:** the `assay`, `fwd`, `rev`, and `probe` columns must **all be
present**. The run fails with `Assay table missing required columns` if any is absent. Only
`probe` may carry an empty *value* (that is how you declare a probe-free assay); `assay`,
`fwd`, and `rev` need real values. The `target_group` and `target_gene` columns are optional
and may be omitted entirely.

**Extra columns:** you may add any additional columns to `assay_table.csv` (e.g.
`reference`, `notes`) to keep your work organized; AssayVal ignores them. The one
exception is `target_gene` — if you include it, it is carried through to
`assay_performance.csv` alongside `fwd`, `rev`, and `probe` to make your life easier.

**Checks run before any search:** assay names must be unique and contain no spaces
(BLAST truncates sequence names at the first space, which would silently orphan the
assay's hits; use underscores), and after modifications are stripped every `fwd`,
`rev` and `probe` must contain only IUPAC letters. Leading/trailing spaces in any cell
are ignored. A failing check stops the run with a message naming the line, assay and
column.

### Editing the assay table in Excel

Excel's plain **"CSV"** format saves in your platform's legacy encoding (Mac Roman on
macOS, Windows-1252 on Windows), and in some European locales uses `;` as the
delimiter. AssayVal reads all of these: a non-UTF-8 file is decoded automatically
(non-ASCII text such as `14–15`, `60 °C` or accented author names survives intact) and
a one-line warning is printed. To avoid the warning, save as **"CSV UTF-8 (Comma
delimited)"**. The `assay_table.csv` written by `assay-val init` starts with a UTF-8
byte-order mark, so Excel opens it as CSV UTF-8 and a plain **Save** keeps that format.
`metadata.csv` is read the same way.

### Oligo notation (LNA and modifications)

Write each oligo the way you would order it. AssayVal interprets every token, keeps
the bases a modification carries, removes the ones that carry no base, and stops with
an error naming any token it does not recognise (rather than guessing, which could
silently delete or add a base).

| Notation | Example | Searched as |
|---|---|---|
| Locked nucleic acid (IDT style) | `TG+GAATC+GTTT+GACTGCATTT` | `TGGAATCGTTTGACTGCATTT`, LNA at positions 3, 8, 12 |
| Locked nucleic acid (brackets) | `TG[G]AATC[G]TTT[G]ACTGCATTT` (or `[+G]`) | same as above |
| IDT 5′ / 3′ modification (no base) | `/56-FAM/ACGT/3IABkFQ/` | `ACGT`. `/5…/` must come first and `/3…/` last. |
| IDT internal quencher / spacer (no base) | `/ZEN/`, `/iTAO/`, `/iSpC3/`, `/iSp9/`, `/iSp18/` | removed |
| Other non-base modification | `ACGT[BHQ1]`, `[AmMC6]` (two or more characters in brackets) | removed |
| 5-methyl-dC | `/iMe-dC/` | `C` |
| Inosine | `/ideoxyI/` or `I` | `N` (scored as matching any base; slightly generous for I·G) |
| Deoxyuridine / uracil | `/ideoxyU/` or `U` | `T` |
| Labelled dT | `/iBiodT/`, `/iFluorT/`, `/iAmMC6T/` | `T` |
| Phosphorothioate bond | `A*C*G*T` | `ACGT` (backbone only) |

`resources/oligos/prep.log` shows, for every oligo, the sequence actually searched,
each modification removed or converted, and LNA positions. The run manifest lists LNA
positions per assay.

**How LNA affects detection calls.** An LNA base pairs like the ordinary base, so
mismatch counting, the 3′-exact rule and the BLAST search are unchanged. But LNA is
used precisely because it sharpens mismatch discrimination: LNA probes are designed so
that a mismatch at the LNA site prevents binding. With `lna_mismatch: exact` (the
default, in `config.yaml`) **a mismatch at any LNA position rejects that primer or probe
hit**, however many mismatches are otherwise allowed. Set `lna_mismatch: count` to
score it as an ordinary mismatch instead. This filter runs after BLAST, so both
settings can be compared on one run with
`assay-val --rescore-from <run> --set lna_mismatch=count`. AssayVal compares sequences
only; it does not model the Tm increase LNA provides, and mismatches *adjacent* to an
LNA (also destabilizing) are counted normally.

**Not supported yet:** RNA bases (`rA`, `[rA]`), 2′-O-methyl RNA (`mA`, `[mA]`) and
rhAmp/rhPCR primers are rejected with an error, because those primers are only
extended after cleavage at the RNA base and need their own model. A lowercase `r`/`m`
is still read as the IUPAC code R/M in an all-lowercase sequence.

> **Excel tip:** a cell that *starts* with `+` (an LNA at the 5′ end, e.g.
> `+TGGAATC…`) is treated by Excel as a formula and becomes `#NAME?`, both when you
> type it and each time Excel re-opens the CSV. Write a 5′-terminal LNA in bracket form
> (`[T]GGAATC…`) instead. AssayVal reports a clear error if it finds `#NAME?` or a
> formula in a sequence cell.

### Per-assay thresholds

The detection thresholds in `config.yaml` apply to every assay. To change them for
individual assays, add any of the columns `max_primer_mismatches`, `prime3_exact_nt`,
`max_probe_mismatches` and `max_amplicon_size` to `assay_table.csv` (`assay-val init`
includes them, blank). A blank cell uses the config value; a whole number overrides
it for that assay only. Typical uses are `max_probe_mismatches: 0` for short,
mismatch-sensitive MGB or LNA probes, and a `max_amplicon_size` that suits each
assay's expected product. Overrides take precedence over `assay-val --set`, so an
assay with its own value is not varied by a threshold sweep. Overrides are listed in
the run manifest and carried into `assay_performance.csv`.

### Probe-free assays

Leave `probe` empty for SYBR/dsDNA-dye chemistry. AssayVal then requires only a valid
amplicon to call a detection: no probe oligo is searched, and the `Primer Only` call is
structurally unreachable for that assay. Any valid amplicon is `Detected`, and
`Not Detected` means no valid amplicon was found. This mirrors the chemistry, where any
double-stranded product fluoresces. Note the consequence for cross-assay comparison,
described under [Interpretation](#interpretation).

Example:

```
assay,probe,fwd,rev,target_group
Assay1,ACGGGACAAAAAGGATGGCGAGTAC,AGCCGAGCGTTACCAGC,CGAACGCAATGATTCTCTGAGC,species:Vibrio mediterranei
Assay2,,GCTACGCCCTCCATCATCC,GCGCGTGATTATCTGATAGC,
```

`Assay1` is probe-based and scored against a target group. `Assay2` has an empty `probe`
(probe-free) and an empty `target_group`, so it is reported but not scored.

---

## Metadata format

`metadata.csv` describes your assemblies. **Only one column is required:**
- `accession` — must match the assembly filenames (`<accession>.fna`).

Every other column is optional and carried through to the per-assembly outputs. To group assemblies in the reports, set `group_by` in `config.yaml` to one or more metadata columns:

```yaml
group_by: "species"                 # single grouping
group_by: ["species", "phenotype"]  # a detection matrix + heatmap per column
```

Grouping-column **names** must not contain whitespace (use e.g. `isolation_source`, not `isolation source`).

If you used `download-assemblies`, the generated `metadata.csv` carries NCBI ANI/BioSample fields and AssayVal **auto-groups by ANI-derived species** when `group_by` is unset. If you bring your own `metadata.csv` without those NCBI columns, you must set `group_by`. When you set `group_by`, assemblies missing from `metadata.csv` or with an empty cell in a grouping column are still analyzed and reported under `Ungrouped`. In the default ANI-auto mode (no `group_by`), only species-confident assemblies get their own group; assemblies that are only genus-resolvable, unclassified, or otherwise low-confidence fold into a single `Unclassified (low confidence ANI)` group. Each assembly's tier is still retained per-assembly as `ani_confidence` (`High`/`Genus`/`Low`) in `detection_by_assembly.csv`, and the run manifest reports the high/genus/low counts.

Because grouping columns are independent, an assay can be scored against a different resolution than the one used to lay out a matrix (e.g. group the report by `species` while scoring a nested clade assay with `target_group: phenotype:protective`).

Example minimal metadata:
```
accession,species,phenotype
GCF_000000001.1,Vibrio mediterranei,protective
GCF_000000002.1,Vibrio harveyi,
```

---

## Detection thresholds

AssayVal reports a detection call per assay per assembly using thresholds set in
`config/config.yaml`. The defaults reflect PCR biochemistry:

| Parameter               | Default | Rationale                                                                                                                                                                                                                         |
| ----------------------- | ------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `max_primer_mismatches` | 2       | A primer with one or two **internal** mismatches still typically primes efficiently; this tolerates strain-level SNPs while excluding poor binders. Counted IUPAC-aware (a degenerate base matches any of its represented bases). |
| `prime3_exact_nt`       | 1       | A mismatch at the 3′-terminal base reliably blocks polymerase extension, so the final base must match exactly regardless of `max_primer_mismatches`. A single 3′-penultimate mismatch usually still amplifies (at reduced efficiency), so it is counted against the mismatch budget instead; raise to 2–3 for stricter extension-critical filtering. |
| `max_probe_mismatches`  | 1       | Used for probe-based assays only. Hydrolysis probes tolerate less mismatch than primers, so the default is stricter.                                                                                                              |
| `max_amplicon_size`     | 500     | Typical qPCR/dPCR amplicons are ~70–200 bp; 500 bp captures valid products while rejecting spurious long-range primer pairings.                                                                                                  |

**Detection calls:**
- `Detected` — valid amplicon found with probe contained within it (probe assays), or valid amplicon found (probe-free assays)
- `Primer Only` — valid amplicon found but probe not detected within it
- `Not Detected` — no valid amplicon found

**Interpreting `Primer Only`:** both primers bind and would amplify, but the probe site is diverged or absent. For a hydrolysis-probe (dPCR/qPCR) assay this usually means **no fluorescent signal** despite amplification, so `assay_performance.csv` counts it as a non-detection; for SYBR/probe-free chemistry any valid amplicon is a detection. `pct_detected_or_primer` in the summaries lets you see both interpretations.

The BLAST search itself is run with deliberately permissive settings
(`evalue=1000`, `perc_identity=70`, `word_size=4`, `max_target_seqs=50000`) so
that no candidate binding site is missed; stringency is enforced downstream by
the mismatch, 3′-exact, and amplicon-size filters above.

Two of those defaults deserve their rationale spelled out, because the obvious
choices are wrong in ways that fail silently:

- **`word_size: 4`, not blastn-short's default of 7.** Every BLAST hit must
  contain at least one exact run of word-length, so an oligo whose mismatches
  break it into runs shorter than the word size is never seeded — it is
  invisible regardless of how permissive `max_primer_mismatches` is. The blind
  spot grows sharply as oligos shorten: at `word_size 7` it hides 7.4% of
  two-mismatch placements in a 17-mer primer, and 7.7% of *single*-mismatch
  placements in a 13-mer MGB probe. `word_size 4` is the minimum blastn allows
  and reduces that to zero for every oligo length and mismatch budget in normal
  use. It is also
  [simulate_PCR](https://doi.org/10.1186/1471-2105-15-237)'s default, for the
  same reason. Raise it only to trade sensitivity for speed.
- **`max_target_seqs: 50000`.** blastn defaults to 500. Because primeval builds
  one database per assembly, the *subjects* are contigs, and a fragmented draft
  assembly can exceed that cap — at which point hits are dropped during the
  search rather than ranked and truncated afterwards
  ([Shah et al. 2019](https://doi.org/10.1093/bioinformatics/bty833)). Setting
  it far above any plausible contig count removes the failure mode.

### Terminal mismatches and BLAST end-trimming

BLAST reports *local* alignments: an alignment is trimmed to the ends that
maximize its score, so it can never extend through a mismatch at an oligo
terminus (with `blastn-short` scoring, +1 match / −3 mismatch, any mismatch
within ~4 bases of either end causes trimming). A naive "100% query coverage"
filter downstream of BLAST would therefore silently discard every hit with a
mismatch near a primer end — no matter how permissive `max_primer_mismatches`
is — producing systematic false negatives for assays whose only template
variation sits at a primer terminus.

AssayVal instead **reconstructs** each gapless hit to full oligo length: the
un-aligned oligo ends are mapped onto the assembly via the hit coordinates,
the complete window is re-read from the genome sequence (reverse-complemented
for minus-strand hits; positions beyond a contig edge are treated as
mismatches), and the IUPAC-aware mismatch count and `prime3_exact_nt` filter
are applied to the full-length window. 5′ mismatches thus count only against
the `max_primer_mismatches` budget — consistent with their minimal effect on
priming — while 3′-terminal mismatches remain governed by the strict
`prime3_exact_nt` rule. Probes are handled the same way against
`max_probe_mismatches`.

---

## Threshold sensitivity analysis (re-scoring)

The [detection thresholds](#detection-thresholds) above are defensible defaults,
not tuned parameters, and a reviewer is entitled to ask whether a reported
sensitivity or specificity is an artifact of where those cutoffs were drawn. The
honest answer is to vary them and show what happens.

Re-running the whole pipeline per threshold setting is the obvious way to do that
and the wrong one: BLAST dominates the cost of a run (roughly 15 MB of transient
output per assembly, ~158 GB for the 10,715-assembly *Vibrionaceae* set), and the
BLAST search does not depend on the thresholds at all — only the scoring does.
`--rescore-from` separates the two. Run BLAST once, keep it, then score it as many
times as you like:

```bash
# 1. one full run, with keep_blast: true set in config.yaml
assay-val --run-name Vpop

# 2. re-score it at a different threshold setting — no BLAST, minutes not hours
assay-val --run-name Vpop_mm1_p3x2 \
  --rescore-from results/Vpop_<date> \
  --set max_primer_mismatches=1 \
  --set prime3_exact_nt=2
```

Each re-score writes a complete, independent run directory — its own
`amplicons/`, `reports/`, heatmaps and `assay_performance.csv` — so the outputs
are directly comparable with the original run. `run_manifest.txt` records the
thresholds used **and** a `rescored_from` line naming the source run, so every
cell of a sweep is self-documenting.

A full grid is a shell loop:

```bash
SRC=results/Vpop_<date>
for mm in 0 1 2 3; do
  for p3 in 1 2 3; do
    assay-val --run-name "Vpop_mm${mm}_p3x${p3}" --rescore-from "$SRC" \
      --set max_primer_mismatches=$mm --set prime3_exact_nt=$p3
  done
done
```

Then collate `results/Vpop_mm*/reports/assay_performance.csv` (each carries its
own `assay`, `sensitivity` and `specificity` columns) to plot sensitivity and
specificity against threshold, one panel per assay.

**Requirements and caveats**

- The source run must have been run with `keep_blast: true`. Without it each
  BLAST TSV is deleted as it is consumed and there is nothing to re-score;
  `--rescore-from` fails immediately and says so.
- **The assemblies must still be present.** Scoring is not a pure function of the
  BLAST output: `detect.py` re-reads each genome to rebuild BLAST end-trimmed
  hits to full oligo length (see [Terminal mismatches and BLAST
  end-trimming](#terminal-mismatches-and-blast-end-trimming)). `assembly_dir`
  must still resolve.
- The work set is taken from the cached BLAST output, not from `assembly_dir`, so
  re-scoring a subsampled run stays subsampled even if the assembly directory has
  since grown.
- `--set` accepts any config key, but changing anything that feeds the BLAST
  search itself (`assay_table`, `blast_*`) is meaningless in re-score mode: those
  results are baked into the cache. Change an assay and you need a full run.
  Per-assay threshold columns and `lna_mismatch` *are* applied at re-score time.
- BLAST output cached by a version before LNA support searched bracket-LNA
  (`[G]`), `/ideoxyI/` and `/iMe-dC/` oligos with those bases deleted. Re-run BLAST
  (a full run) for any assay written that way rather than re-scoring its old cache.
- Re-scoring is not free — the detection step still runs per assembly — but it
  skips the database build and the search, which is where the time and disk go.

---

## Outputs

All outputs are written to the run directory (`results/<run-name>_<date>/`):

```
results/<run-name>_<date>/
├── amplicons/
│   └── {accession}.csv           # per-assembly, one row per assay: detection call,
│                                  #   mismatch counts, amplicon sizes/contigs/positions,
│                                  #   and (if enabled) amplicon sequences
├── amplicon_fasta/
│   ├── {assay}_amplicons.fasta   # one per assay: every amplicon (+ flanks) — see below
│   └── amplicon_index.csv        # one row per FASTA record
├── run.log                       # full run log
└── reports/
    ├── {column}_detection_matrix.csv   # one per group_by column: group × assay % detected
    ├── detection_summary_long.csv       # tidy long table: assay × grouping × group (incl. misses)
    ├── assay_summary.xlsx               # same as detection_summary_long, one worksheet per assay
    ├── assay_performance.csv            # per-assay sensitivity / specificity vs target_group
    ├── detection_by_assembly.csv        # one row per assembly: metadata + 0/1 per assay
    ├── run_manifest.txt                 # parameters, tool versions, checksums, input accounting
    ├── per_assay/                       # {assay}_results.csv — every assembly's call for one assay
    └── figures/
        └── {column}_detection_heatmap.pdf/png   # one per group_by column
```

**Notes**
- `{column}_detection_matrix.csv` and its matching `figures/{column}_detection_heatmap.*`
  are written once per `group_by` column (or once as `group_detection_matrix.csv` /
  `group_detection_heatmap.*` in ANI-auto mode, since the auto-derived column is named
  `group`). See [Metadata format](#metadata-format).
- `detection_summary_long.csv` (one tidy table) and `assay_summary.xlsx` (one worksheet
  per assay) are the **same data in different layouts**, not a duplicate file.
- `detection_summary_long.csv` includes assay × grouping-column × group combinations an
  assay *misses* (`pct_detected = 0`), so you can see both what an assay detects and
  what it doesn't, across every `group_by` column.
- `assay_performance.csv` scores each assay against its declared `target_group`. Because
  targets are `column:value`, a nested clade assay (`phenotype:protective`) and a
  species-wide assay (`species:Vibrio mediterranei`) are each scored against the right
  axis in one run. See [Interpreting `assay_performance.csv`](#interpreting-assay_performancecsv)
  for the full column reference and how to read the numbers.
- `n_multi_amplicon` / `max_amplicons` flag assemblies where an assay produces more than
  one product. Interpretation is target-dependent: for a single-copy target this can
  indicate off-target priming or overestimated abundance, whereas for a multi-copy
  target (e.g. 16S rRNA) multiple amplicons per genome are expected.
- `detection_by_assembly.csv` joins the binary detection matrix (`1` = Detected) with
  every column of the input `metadata.csv`, one row per assembly.
- Raw BLAST output and the per-assembly `logs/` directory are intermediate and are
  deleted once a run completes successfully. Set `keep_blast: true` / `keep_logs: true`
  in `config.yaml` to retain them; the top-level `run.log` is always kept.

### Amplicon FASTA files

`amplicon_fasta/` holds one `{assay}_amplicons.fasta` per assay containing **every
valid amplicon** that assay produced — `Detected` and `Primer Only` alike — ready to
import into Geneious or an aligner, e.g. to design a probe for a probe-free assay, or to
see which primer/probe positions vary in off-target taxa. Every assay gets a file; an
assay with no amplicons gets an empty one.

- **Orientation:** every record reads from the forward primer to the reverse primer
  (amplicons on the minus strand are reverse-complemented), so records align directly.
- **Case:** the primer-to-primer amplicon is UPPERCASE; `amplicon_flank_bp` (default 50)
  of flanking sequence on each side is lowercase. Flanks are clipped at contig ends.
  Primer and probe sites are the genome's own sequence, so mismatches are visible.
  (Some aligners, e.g. MAFFT, lowercase everything unless told to preserve case.)
- **Names:** `{accession}_amp{k}of{n}_{organism}`, e.g.
  `GCF_000295695.2_amp3of10_Bacillus_anthracis_str._BF1`. Every amplicon is numbered,
  so assemblies with several copies (e.g. 16S rRNA) stay unique, and `k` follows the
  order of `amplicon_starts` / `contig_ids` in `amplicons/{accession}.csv`. Geneious
  uses this as the sequence name, so alignment rows show taxa; the accession and
  organism remain searchable. The organism is `organism_name` from `metadata.csv`
  (falling back to the first `group_by` column, then `unknown`).
- **Description:** the rest of the header line, e.g.
  `organism="Bacillus anthracis str. BF1" assay=Bac_16S_Food assembly_call=Detected
  probe=yes contig=NZ_CP047131.1:82367-82426(+) amplicon_bp=60 flank_bp=50,50 fwd_mm=0
  rev_mm=0 probe_mm=0`. `probe` is per copy (`yes`/`no`; `none` for probe-free
  assays): in a `Detected` assembly, some copies may still lack the probe site.
- `amplicon_index.csv` lists every record (its FASTA name, file and all header fields)
  in one table.

Set `amplicon_fasta: false` in `config.yaml` to skip the export.

---

## Interpreting `assay_performance.csv`

For every assay carrying a `target_group`, AssayVal classifies **each assembly in the run**
two ways and cross-tabulates them:

- **Truth** — is this assembly in the assay's intended target group? (i.e. does
  `metadata[target_column]` equal the value in `target_group`?)
- **Call** — did the assay produce a detection in this assembly? (`detection_call` is `Detected`)

|                         | Detected (in silico)  | Not detected           |
| ----------------------- | --------------------- | ---------------------- |
| **In target group**     | `tp` (correct hit)    | `fn` (missed target)   |
| **Not in target group** | `fp` (off-target hit) | `tn` (correct reject)  |

### Columns

| Column             | Meaning                                                                                                                            |
| ------------------ | ---------------------------------------------------------------------------------------------------------------------------------- |
| `assay`            | Assay name, from `assay_table.csv`.                                                                                                |
| `fwd`, `rev`, `probe` | Carried over verbatim from `assay_table.csv` (including any `/ZEN/`-style modification notation), so the results file is self-contained. |
| `target_gene`      | Carried over from `assay_table.csv` when you include that optional column.                                                          |
| `target_group`     | The target exactly as written in the assay table (e.g. `phenotype:protective`).                                                    |
| `target_column`    | The metadata column the target resolved to (e.g. `phenotype`). A bare value resolves to the primary `group_by` column.             |
| `n_target`         | Assemblies in the target group (`tp + fn`).                                                                                        |
| `n_nontarget`      | Assemblies outside the target group (`fp + tn`).                                                                                   |
| `tp`               | **True positives** — target assemblies the assay detected.                                                                          |
| `fn`               | **False negatives** — target assemblies the assay failed to detect (i.e. sensitivity gaps).                                        |
| `fp`               | **False positives** — non-target assemblies the assay detected (i.e. cross-reactivity).                                            |
| `tn`               | **True negatives** — non-target assemblies the assay correctly did not detect.                                                     |
| `sensitivity`      | `100 × tp / (tp + fn)` — the % of intended targets the assay detects (assay **inclusivity**).                                       |
| `specificity`      | `100 × tn / (tn + fp)` — the % of non-targets the assay correctly rejects (assay **exclusivity**).                                  |
| `n_detected_total` | Every assembly this assay detected in the run (`tp + fp`).                                                                          |

Control/reference assays with a blank `target_group` still get a row, with blank metrics and
only `n_detected_total` filled, so the file remains a complete roster of every assay. Any rate
whose denominator is zero is left blank rather than reported as `0` — e.g. a `target_group`
that matches no assemblies yields `n_target = 0`, a blank `sensitivity`, and a warning in
`run.log`.

### Interpretation

- **These are in silico predictions, not wet-lab performance.** They report whether the primers
  and probe have acceptable binding sites in each genome under the configured
  [detection thresholds](#detection-thresholds). AssayVal does NOT test amplification efficiency, Tm, secondary
  structure, or partitioning behaviour. Read them as the sequence-level expectation that
  empirical validation is tested against.
- **Every denominator is your input assembly set.** `specificity = 100` means "no off-target
  detections *among the genomes you supplied*," not "no off-targets exist." See
  [Assay specificity validation](#assay-specificity-validation).
- **`Primer Only` counts as a non-detection** (see [Detection thresholds](#detection-thresholds)),
  so a probe-based assay that amplifies off-target *without* probe binding is correctly **not**
  counted as a false positive.
- **Probe-based and probe-free assays are not directly comparable.** A probe-free
  ([SYBR/dsDNA-dye](#probe-free-assays)) assay scores a detection on two binding sites; a
  probe-based assay needs three. All else equal that gives probe-free assays systematically
  *higher* apparent sensitivity and *lower* apparent specificity. The `probe` column carried into this file tells you
  which is which (blank = probe-free), so compare like with like.

---

## Downloading assemblies

`download-assemblies` is installed with the environment and runs from any
directory.

### Basic usage

```bash
download-assemblies -t "Taxon name" -o assemblies/
```

Options:

| Flag | Description | Default |
|------|-------------|---------|
| `-t TAXON` | Taxon name or NCBI tax ID (required; **repeatable** — see below) | — |
| `-o OUTDIR` | Output directory | `assemblies` |
| `-l LEVELS` | Assembly levels (comma-separated) | `complete,chromosome,scaffold,contig` |
| `-s SOURCE` | Assembly source: `refseq`, `genbank`, or `all` | `refseq` |
| `-n N` | Randomly keep at most N assemblies per `-t` taxon | all |
| `-r SEED` | Integer random seed for `-n` selection | `0` |
| `-e EMAIL` | NCBI e-mail (or set `NCBI_EMAIL` env var) | — |
| `-k API_KEY` | NCBI API key for higher rate limits (or set `NCBI_API_KEY` env var) | — |

### Multiple taxa

Pass `-t` more than once to download the **de-duplicated union** of several taxa
into a single output directory:

```bash
download-assemblies \
    -t "Vibrio jasicida" \
    -t "Vibrio owensii" \
    -o assemblies/
```

Each taxon is queried in turn; assemblies shared between taxa (e.g. when one
query nests inside another) are downloaded and written to `metadata.csv` only
once. Mixing names and tax IDs is fine, e.g. `-t "Vibrio owensii" -t 661487`.

### Setting your NCBI API key once

An [NCBI API key](https://www.ncbi.nlm.nih.gov/account/) raises your download
rate limit from 3 to 10 requests/sec. This is worth setting if you're downloading large sets of assemblies (i.e., hundreds or thousands of assemblies). Rather than passing `-k` every time, save it once:

```bash
cp config/ncbi_credentials.example.sh config/ncbi_credentials.sh
# edit config/ncbi_credentials.sh and paste your key into NCBI_API_KEY
```

`download-assemblies` sources this file automatically on every run.
To keep the key elsewhere, point `ASSAYVAL_CREDENTIALS` at your own file.

The key is resolved as: **`-k` flag → `NCBI_API_KEY` environment variable →
credentials file** (first one set wins).

### Scaling down large taxa

By default `download-assemblies` fetches every assembly matching your level and
source filters. To cap how many are downloaded **per `-t` taxon**, add `-n N`
(and, optionally, `-r SEED` — an integer seed, default `0`, for reproducible
selection):

```bash
download-assemblies -t Klebsiella -n 20 -r 1 -o assemblies/
```

Because the downloader is resume-aware and de-duplicates into one output
directory, mix different per-taxon counts by composing invocations into the same
`-o` directory:

```bash
download-assemblies -t Klebsiella          -n 20 -r 1 -o assemblies/
download-assemblies -t "Escherichia coli"  -n 20 -r 1 -o assemblies/
download-assemblies -t Vibrionaceae               -o assemblies/   # all
# -> assemblies/ holds 20 Klebsiella + 20 E. coli + all Vibrionaceae
```

Selection is uniform-random within your `-l`/`-s` filters. The seed makes it
deterministic *for a given NCBI result set*; because NCBI's holdings grow over
time, the authoritative record of what was actually downloaded is the run's
`metadata.csv` and `download.log` (which reports `sampled N of M available` per
taxon). `-n` never errors when a taxon has fewer than `N` assemblies — it simply
takes all of them. `-n` is applied to each `-t` taxon **before** cross-taxon
de-duplication, so with nested or overlapping taxa the combined set can end up
smaller than N × (number of taxa).

Composing invocations into one `-o` directory accumulates both the `.fna` files
**and** `metadata.csv` — each run merges its metadata into the existing CSV (by
accession) rather than overwriting it, so the recipe above's final `metadata.csv`
covers every assembly downloaded across all three invocations.

**Curated sets:** subsampling is for scaling down broad taxa, not for building a
specific positive set. If you need particular genomes (e.g. an allele-diverse set
of a resistance gene's carriers), download them directly and drop the `.fna`
files into your `assembly_dir` — AssayVal reads any `.fna`:

```bash
datasets download genome accession GCF_XXXXXXXXX.1 --include genome
# unzip and place the .fna in assemblies/
```

### HPC / SLURM

Local available storage requirements for AssayVal are directly scaled by the assembly dataset provided (i.e., you need space to store the downloaded assemblies you provide AssayVal!). If needed, AssayVal runs can easily be submitted in a SLURM environment using the wrapper below.

Wrap the command in an sbatch job for large downloads:

```bash
sbatch --time=24:00:00 --mem=8G \
  --wrap="download-assemblies -t Vibrionaceae -o assemblies/"
```

`download-assemblies` is resume-aware: if interrupted, re-running it will skip assemblies already successfully downloaded.

---

## Assay specificity validation

AssayVal tests assay sensitivity and specificity against **your input assembly dataset**. The scope of specificity testing is therefore determined by which assemblies you provide.

**Recommended workflow for specificity screening:**

1. **Single-primer BLAST screen** (NCBI web interface): Individually BLAST each primer and probe sequence against the NCBI `nt` database, *excluding* your target taxon. This identifies any off-target binding sites outside your group of interest. If no hits are returned for any oligo, off-target amplification outside the taxon is extremely unlikely (a primer must bind for any amplicon to form).

2. **Expand the input dataset if needed**: If step 1 returns hits in a non-target taxon, download assemblies from that taxon and add them to your `assembly_dir`. AssayVal will then determine whether those single-primer hits form complete, detectable amplicons.

This two-stage approach is computationally efficient, such that you only download and evaluate assemblies in taxa where off-target primer binding is possible.

Step 1 is a manual pre-screen performed through the NCBI web interface; it is not part of the reproducible AssayVal pipeline. The reproducible specificity evidence for a manuscript comes from AssayVal's `assay_performance.csv` over your assembled dataset (Step 2 onward).

> **In the future:** A future release will support BLASTing directly against NCBI pre-built reference databases (e.g., `ref_prok_rep_genomes`) as a single-step broader specificity check, without requiring manual assembly downloads.

---

## Test dataset

A small set of 5 assemblies is provided for validating your installation:

```bash
bash test_data/download_test_data.sh
```

Then update `config/config.yaml`:

```yaml
assembly_dir: "test_data/assemblies"
metadata: "test_data/assemblies/metadata.csv"
```

And run from the repo root (using the repo's config directly):

```bash
assay-val --run-name test --configfile config/config.yaml
```

The test dataset covers all detection scenarios: `Detected` (including via minus-strand primer binding), `Primer Only`, and `Not Detected`.

### Example output

With `group_by: "species"` set and a `target_group` of `species:Vibrio mediterranei`
declared for assay `VmedA`, the report format looks like this:

```
# species_detection_matrix.csv
group,n_assemblies,VhPath,VmedA
Vibrio harveyi,2,100.0,0.0
Vibrio mediterranei,3,0.0,66.67

# assay_performance.csv
assay,fwd,rev,probe,target_group,target_column,n_target,n_nontarget,tp,fn,fp,tn,sensitivity,specificity,n_detected_total
VmedA,GCTACGCCC…,GCGCGTGAT…,ACGACCTTC…,species:Vibrio mediterranei,species,3,2,2,1,0,2,66.67,100.0,2
```

---

## assay-design — clade-specific target discovery

**assay-design** is the companion tool used to *design* the assays that AssayVal
validates. It parses a [Panaroo](https://gtonkinhill.github.io/panaroo/)
pangenome to find orthologs that are conserved within a clade and specific to it
(absent elsewhere), then extracts representative protein/nucleotide sequences —
the candidate targets from which the Vpop dPCR assays were built.

It installs its own command, `assay-design`, on your PATH:

```bash
assay-design \
    --matrix          gene_presence_absence.csv \
    --isolates-dir    isolate_groups \
    --gene-data       gene_data.csv \
    --representatives representatives.tsv \
    --output-dir      results
```

A tiny, worked example lives in
[`assay-design/example/`](assay-design/example/). For inputs, thresholds, and
full usage, see the [assay-design README](assay-design/README.md).

---

## Citation

If you use AssayVal, please cite the archived release (TBD!):

> Smith S, et al. (2026) *[manuscript title]*. *[journal]*. doi:[doi]
> AssayVal [version] (2026). Zenodo. doi:[zenodo-doi]

---

## License

MIT — see [LICENSE](LICENSE).
