# Download Subsampling Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a per-taxon random cap (`-n N`, `-r SEED`) to `download-assemblies` so users can scale down large taxa, without changing the default "download all" behavior.

**Architecture:** A new stdlib-only Python helper `scripts/download/subsample.py` holds the testable selection logic (`sample_records`); `download_assemblies.sh` samples each `-t` taxon's `datasets summary` output through it inside the existing query loop, before the unchanged cross-taxon de-duplication and download steps.

**Tech Stack:** Bash (macOS bash 3.2-compatible), Python 3 standard library, pytest.

## Global Constraints

- **Do NOT commit.** The user commits manually. End each task by running the test suite and leaving changes in the working tree — replace any "commit" step with "run the suite, confirm green."
- Run tests with the project env on PATH: `PATH="$HOME/miniforge3/envs/assayval/bin:$PATH" python -m pytest -q`. All **128 existing tests must stay green** at every task boundary.
- `scripts/download/subsample.py` uses the **Python standard library only** (match `scripts/download/parse_metadata.py`).
- Preserve **macOS bash 3.2** compatibility in `download_assemblies.sh` (no `mapfile`; guard `set -u` array expansion as the existing code does; `[[ … =~ … ]]` is fine in 3.2).
- New CLI flags are **short-flag `getopts`** style matching the existing `-t/-o/-l/-s/-e/-k`: `-n N` (per-taxon cap) and `-r SEED` (default `0`).
- **Determinism contract:** `sample_records` sorts by `accession` then applies a seeded shuffle, so the selection is independent of the order NCBI returned records in.
- No `-n` → behavior is byte-for-byte the current behavior (sampling code is inert).

---

## Task 1: `subsample.py` — testable selector + thin CLI

**Files:**
- Create: `scripts/download/subsample.py`
- Create test: `tests/test_subsample.py`

**Interfaces:**
- Produces: `sample_records(records: list[dict], n: int, seed: int) -> list[dict]` — returns up to `n` records; sorts by `accession`, seeded-shuffles, truncates; `n <= 0` or `n >= len(records)` returns all (sorted); every returned record is an input record.
- Produces CLI: `python3 subsample.py --jsonl IN --n N --seed S --out OUT` — reads JSON-Lines, writes the sampled subset as JSON-Lines.

- [ ] **Step 1: Write the failing tests** — create `tests/test_subsample.py`:

```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts" / "download"))

from subsample import sample_records


def _recs(*accs):
    return [{"accession": a, "extra": i} for i, a in enumerate(accs)]


def test_deterministic_same_seed():
    recs = _recs(*[f"GCF_{i:03d}" for i in range(50)])
    a = [r["accession"] for r in sample_records(recs, 10, seed=1)]
    b = [r["accession"] for r in sample_records(recs, 10, seed=1)]
    assert a == b and len(a) == 10


def test_order_independent():
    import random as _r
    recs = _recs(*[f"GCF_{i:03d}" for i in range(50)])
    shuffled = recs[:]; _r.Random(99).shuffle(shuffled)
    a = [r["accession"] for r in sample_records(recs, 10, seed=1)]
    b = [r["accession"] for r in sample_records(shuffled, 10, seed=1)]
    assert a == b


def test_seed_sensitivity():
    recs = _recs(*[f"GCF_{i:03d}" for i in range(200)])
    a = set(r["accession"] for r in sample_records(recs, 20, seed=1))
    b = set(r["accession"] for r in sample_records(recs, 20, seed=2))
    assert a != b


def test_n_ge_len_returns_all():
    recs = _recs("GCF_002", "GCF_001", "GCF_003")
    got = sample_records(recs, 10, seed=0)
    assert set(r["accession"] for r in got) == {"GCF_001", "GCF_002", "GCF_003"}


def test_n_nonpositive_returns_all():
    recs = _recs("GCF_001", "GCF_002")
    assert len(sample_records(recs, 0, seed=0)) == 2
    assert len(sample_records(recs, -5, seed=0)) == 2


def test_subset_and_size():
    recs = _recs(*[f"GCF_{i:03d}" for i in range(30)])
    got = sample_records(recs, 7, seed=3)
    assert len(got) == 7
    inputs = set(r["accession"] for r in recs)
    assert all(r["accession"] in inputs for r in got)


def test_cli_roundtrip(tmp_path):
    import subprocess, sys as _sys, json
    inp = tmp_path / "in.jsonl"
    inp.write_text("\n".join(json.dumps({"accession": f"GCF_{i:03d}"}) for i in range(40)) + "\n")
    out = tmp_path / "out.jsonl"
    script = Path(__file__).parent.parent / "scripts" / "download" / "subsample.py"
    subprocess.run([_sys.executable, str(script), "--jsonl", str(inp),
                    "--n", "5", "--seed", "1", "--out", str(out)], check=True)
    lines = [l for l in out.read_text().splitlines() if l.strip()]
    accs = [json.loads(l)["accession"] for l in lines]
    assert len(lines) == 5 and len(set(accs)) == 5
```

- [ ] **Step 2: Run to confirm failure**

Run: `PATH="$HOME/miniforge3/envs/assayval/bin:$PATH" python -m pytest tests/test_subsample.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'subsample'`.

- [ ] **Step 3: Implement** — create `scripts/download/subsample.py`:

```python
#!/usr/bin/env python3
"""subsample.py

Randomly select up to N assembly-summary records, deterministically for a given
seed. Used by download_assemblies.sh to cap the number of assemblies kept per
taxon. Reads NCBI `datasets summary genome ... --as-json-lines` output (one JSON
record per line) and writes the sampled subset as JSON-Lines.

Usage:
    python3 subsample.py --jsonl in.jsonl --n 20 --seed 0 --out sampled.jsonl
"""
import argparse
import json
import random
import sys


def sample_records(records, n, seed):
    """Return up to n records, chosen deterministically for a given seed.

    Records are sorted by 'accession' (canonical order), then shuffled with a
    seeded RNG and truncated to n. n <= 0 or n >= len(records) returns all
    records. The sort makes the result independent of input order. Every
    returned record is one of the inputs.
    """
    ordered = sorted(records, key=lambda r: r.get("accession", ""))
    if n <= 0 or n >= len(ordered):
        return ordered
    shuffled = ordered[:]
    random.Random(seed).shuffle(shuffled)
    return shuffled[:n]


def _read_jsonl(path):
    records = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    return records


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jsonl", required=True, help="input JSON-Lines summary")
    ap.add_argument("--n", type=int, required=True, help="max records to keep")
    ap.add_argument("--seed", type=int, default=0, help="random seed (default: 0)")
    ap.add_argument("--out", required=True, help="output JSON-Lines path")
    args = ap.parse_args(argv)

    records = _read_jsonl(args.jsonl)
    kept = sample_records(records, args.n, args.seed)
    with open(args.out, "w") as out:
        for r in kept:
            out.write(json.dumps(r) + "\n")
    print(f"subsample: kept {len(kept)} of {len(records)} (n={args.n}, seed={args.seed})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
```

- [ ] **Step 4: Run to confirm pass**

Run: `PATH="$HOME/miniforge3/envs/assayval/bin:$PATH" python -m pytest tests/test_subsample.py -q`
Expected: PASS (7 passed).

- [ ] **Step 5: Gate** — full suite green: `PATH="$HOME/miniforge3/envs/assayval/bin:$PATH" python -m pytest -q` → `135 passed` (128 + 7). Leave changes unstaged (do not commit).

---

## Task 2: Wire `-n`/`-r` into `download_assemblies.sh`

**Files:**
- Modify: `scripts/download/download_assemblies.sh`
- Test: `tests/test_download_cli.py`

**Interfaces:**
- Consumes: `scripts/download/subsample.py` (Task 1), invoked as `python3 "${SUBSAMPLE}" --jsonl … --n … --seed … --out …`.
- Produces: `download-assemblies -n N -r SEED` — per-taxon random cap, validated, logged.

- [ ] **Step 1: Write the failing tests** — add to `tests/test_download_cli.py`:

```python
def test_help_lists_subsample_flags():
    exe = shutil.which("download-assemblies")
    assert exe
    out = subprocess.run([exe, "-h"], capture_output=True, text=True)
    combined = out.stdout + out.stderr
    assert "-n" in combined and "-r" in combined
    assert "seed" in combined.lower()


def test_invalid_n_rejected(tmp_path):
    exe = shutil.which("download-assemblies")
    assert exe
    r = subprocess.run([exe, "-t", "Foo", "-n", "abc", "-o", str(tmp_path / "asm")],
                       capture_output=True, text=True)
    assert r.returncode != 0
    assert "-n" in (r.stdout + r.stderr)


def test_n_zero_rejected(tmp_path):
    exe = shutil.which("download-assemblies")
    assert exe
    r = subprocess.run([exe, "-t", "Foo", "-n", "0", "-o", str(tmp_path / "asm")],
                       capture_output=True, text=True)
    assert r.returncode != 0
```

- [ ] **Step 2: Run to confirm failure**

Run: `PATH="$HOME/miniforge3/envs/assayval/bin:$PATH" python -m pytest tests/test_download_cli.py -q`
Expected: FAIL (help lacks `-n`/`-r`; bad `-n` currently ignored so exit is 0 or reaches network).

- [ ] **Step 3a: Document the flags** — in `download_assemblies.sh`, in the `# Options:` header comment block, add after the `-s SOURCE` line:

```bash
#   -n N           Randomly keep at most N assemblies per -t taxon (default: all)
#   -r SEED        Integer random seed for -n selection (default: 0)
```

- [ ] **Step 3b: Add defaults** — after the `ASSEMBLY_SOURCE="refseq"` default line:

```bash
LIMIT=""       # empty = keep all (no subsampling)
SEED=0
```

- [ ] **Step 3c: Parse the flags** — change the getopts string and add cases:

```bash
while getopts "t:o:l:s:n:r:e:k:h" opt; do
    case $opt in
        t) TAXA+=("$OPTARG") ;;
        o) OUTDIR="$OPTARG" ;;
        l) LEVELS="$OPTARG" ;;
        s) ASSEMBLY_SOURCE="$OPTARG" ;;
        n) LIMIT="$OPTARG" ;;
        r) SEED="$OPTARG" ;;
        e) EMAIL="$OPTARG" ;;
        k) API_KEY="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done
```

- [ ] **Step 3d: Validate** — immediately after the existing `-s` validation `case … esac` block, add:

```bash
if [[ -n "${LIMIT}" ]]; then
    if ! [[ "${LIMIT}" =~ ^[0-9]+$ ]] || [[ "${LIMIT}" -lt 1 ]]; then
        echo "ERROR: -n must be a positive integer (got '${LIMIT}')." >&2
        exit 1
    fi
fi
if ! [[ "${SEED}" =~ ^-?[0-9]+$ ]]; then
    echo "ERROR: -r (seed) must be an integer (got '${SEED}')." >&2
    exit 1
fi
```

- [ ] **Step 3e: Sample per taxon in the query loop** — replace the whole block:

```bash
# Query each taxon and accumulate the raw summaries. Overlapping taxa are
# de-duplicated by accession in the next step.
: > "${SUMMARY_ALL}"
for TAX in "${TAXA[@]}"; do
    log "  querying: ${TAX}"
    datasets summary genome taxon "${TAX}" \
        --assembly-source "${ASSEMBLY_SOURCE}" \
        --assembly-level "${LEVELS}" \
        --as-json-lines \
        ${APIKEY_ARGS[@]+"${APIKEY_ARGS[@]}"} >> "${SUMMARY_ALL}"
done
```

with:

```bash
# Query each taxon. When -n is set, randomly cap each taxon's summaries via
# subsample.py before accumulating. Overlapping taxa are de-duplicated by
# accession in the next step (so a genome shared between taxa is kept once).
SUBSAMPLE="${SCRIPT_DIR}/subsample.py"
: > "${SUMMARY_ALL}"
for TAX in "${TAXA[@]}"; do
    log "  querying: ${TAX}"
    TAX_JSONL="${WORK_TMP}/taxon_summary.jsonl"
    datasets summary genome taxon "${TAX}" \
        --assembly-source "${ASSEMBLY_SOURCE}" \
        --assembly-level "${LEVELS}" \
        --as-json-lines \
        ${APIKEY_ARGS[@]+"${APIKEY_ARGS[@]}"} > "${TAX_JSONL}"
    AVAIL=$(grep -c . "${TAX_JSONL}" 2>/dev/null || echo 0)
    if [[ -n "${LIMIT}" ]]; then
        TAX_SAMPLED="${WORK_TMP}/taxon_sampled.jsonl"
        python3 "${SUBSAMPLE}" --jsonl "${TAX_JSONL}" --n "${LIMIT}" \
            --seed "${SEED}" --out "${TAX_SAMPLED}"
        KEPT=$(grep -c . "${TAX_SAMPLED}" 2>/dev/null || echo 0)
        log "  ${TAX}: sampled ${KEPT} of ${AVAIL} available (seed=${SEED})"
        cat "${TAX_SAMPLED}" >> "${SUMMARY_ALL}"
    else
        log "  ${TAX}: ${AVAIL} available (no subsampling)"
        cat "${TAX_JSONL}" >> "${SUMMARY_ALL}"
    fi
done
```

(`SCRIPT_DIR` and `WORK_TMP` are already defined earlier in the script; `subsample.py` lives beside `parse_metadata.py`, which `SCRIPT_DIR` already points at.)

- [ ] **Step 3f: Log the settings** — after the existing `log "Assembly source: ${ASSEMBLY_SOURCE}"` line, add:

```bash
if [[ -n "${LIMIT}" ]]; then
    log "Subsampling  : up to ${LIMIT} per taxon (seed=${SEED})"
fi
```

- [ ] **Step 4: Run to confirm pass**

Run: `PATH="$HOME/miniforge3/envs/assayval/bin:$PATH" python -m pytest tests/test_download_cli.py -q`
Expected: PASS (all, including the three new tests).

- [ ] **Step 5: Gate** — full suite green: `PATH="$HOME/miniforge3/envs/assayval/bin:$PATH" python -m pytest -q` → `138 passed` (135 + 3). Leave changes unstaged.

---

## Task 3: Document "Scaling down large taxa" in the README

**Files:**
- Modify: `README.md`

Gate: self-review that the docs match the flags implemented in Task 2. No tests.

- [ ] **Step 1: Add the subsection.** In the `## Downloading assemblies` section, immediately **before** the `### HPC / SLURM` subsection, insert:

````markdown
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
takes all of them.

**Curated sets:** subsampling is for scaling down broad taxa, not for building a
specific positive set. If you need particular genomes (e.g. an allele-diverse set
of a resistance gene's carriers), download them directly and drop the `.fna`
files into your `assembly_dir` — AssayVal reads any `.fna`:

```bash
datasets download genome accession GCF_XXXXXXXXX.1 --include genome
# unzip and place the .fna in assemblies/
```
````

- [ ] **Step 2: Gate** — reread the subsection against `download_assemblies.sh`: confirm `-n`, `-r`, the default seed `0`, and the "takes all when fewer than N" behavior all match the code. Fix any drift. Full suite still `138 passed` (docs don't affect tests).

---

## Self-Review

- **Spec coverage:** `-n`/`-r` flags + per-taxon semantics → Task 2; `sample_records` sort-then-seeded-shuffle determinism → Task 1; validation (`n<=0`/non-integer) → Task 2 Step 3d + tests; `N >= available` → all → Task 1 `test_n_ge_len_returns_all`; per-taxon sampling in loop + global dedup unchanged → Task 2 Step 3e; logging seed + sampled/available → Task 2 Steps 3e/3f; reproducibility note → Task 3; unit tests (determinism, seed sensitivity, subset, size, n edge cases) → Task 1; CLI `-h`/bad-`-n` tests → Task 2; docs incl. side-load pointer → Task 3. All spec sections covered.
- **Placeholder scan:** none — every code and test block is complete.
- **Type consistency:** `sample_records(records, n, seed)` is defined in Task 1 and called identically by the Task 2 CLI (`--jsonl/--n/--seed/--out`). `LIMIT`/`SEED`/`SUBSAMPLE`/`TAX_JSONL`/`TAX_SAMPLED` shell variables are consistent across Task 2 steps.

---

## Execution Handoff

Plan saved to `docs/superpowers/plans/2026-07-28-download-subsampling.md`.
