#!/usr/bin/env python3
"""amrfinder_status.py

OPTIONAL AssayVal companion helper — NOT part of the pipeline or the conda
environment, and not required by AssayVal.

Run NCBI AMRFinderPlus over a directory of genome assemblies and append
gene-presence "status" columns (e.g. kpc_status, oxa48_status) to a AssayVal
metadata.csv. Those columns become the ground truth for scoring an AMR-gene
assay in AssayVal via `target_group = <column>:positive` (see the AssayVal
README, "Interpreting assay_performance.csv").

Requires AMRFinderPlus installed separately (`conda install -c bioconda
ncbi-amrfinderplus`) with its database set up (`amrfinder -u`). AMRFinderPlus
is intentionally NOT a AssayVal dependency — AssayVal only cares about the
metadata column, not how you populate it.

Each --gene defines one status column as NAME=PATTERN:

    --gene kpc_status=blaKPC
    --gene "oxa48_status=blaOXA-(48|181|232|204|244)"

NAME is the new metadata column. PATTERN is a case-insensitive regular
expression matched against AMRFinderPlus 'Gene symbol' values for each assembly.
The column is set to "positive" if any detected gene symbol matches, else
"negative". Assemblies with no matching .fna, or whose AMRFinderPlus run failed,
are left blank (unscored) so they are not silently counted as negatives.

NOTE on OXA-48: AMRFinderPlus reports the specific allele (blaOXA-48,
blaOXA-181, blaOXA-232, ...). Decide whether your assay targets ONLY blaOXA-48
(use an anchored pattern like `blaOXA-48$`) or the OXA-48-LIKE family (list the
alleles your assay is meant to detect). This choice defines what counts as a
true positive vs a cross-reactivity false positive.

AMRFinderPlus is run WITHOUT --organism: acquired carbapenemase genes
(blaKPC / blaOXA-48-like) are detected without it, and a single organism can't
be set for an assembly set spanning many genera. (--organism only matters for
chromosomal point-mutation resistance, which this helper does not target.)

Usage:
    python3 amrfinder_status.py \\
        --assembly-dir assemblies/ \\
        --metadata metadata.csv \\
        --gene kpc_status=blaKPC \\
        --gene "oxa48_status=blaOXA-(48|181|232|204|244)"
"""
import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path


def parse_gene_specs(specs):
    """['name=regex', ...] -> [(name, compiled_regex), ...]."""
    out = []
    for s in specs:
        if "=" not in s:
            raise ValueError(f"--gene must be NAME=PATTERN (got '{s}')")
        name, pattern = s.split("=", 1)
        name = name.strip()
        if not name:
            raise ValueError(f"--gene has an empty column name (got '{s}')")
        try:
            rx = re.compile(pattern, re.IGNORECASE)
        except re.error as e:
            raise ValueError(f"--gene '{name}' has an invalid regex '{pattern}': {e}")
        out.append((name, rx))
    return out


def gene_symbols_from_tsv(tsv_path):
    """Return the set of 'Gene symbol' values from an AMRFinderPlus TSV.

    Robust to the column-name variation across AMRFinderPlus versions
    ('Gene symbol' historically, 'Element symbol' in newer releases).
    """
    symbols = set()
    with open(tsv_path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header:
            return symbols
        idx = None
        for col in ("Gene symbol", "Element symbol", "Gene"):
            if col in header:
                idx = header.index(col)
                break
        if idx is None:
            return symbols
        for row in reader:
            if len(row) > idx and row[idx]:
                symbols.add(row[idx])
    return symbols


def status_for(symbols, gene_specs):
    """{column_name: 'positive'|'negative'} for one assembly's gene symbols."""
    return {
        name: ("positive" if any(rx.search(sym) for sym in symbols) else "negative")
        for name, rx in gene_specs
    }


def amrfinder_version():
    try:
        r = subprocess.run(["amrfinder", "--version"], check=True,
                           capture_output=True, text=True)
        return (r.stdout.strip() or r.stderr.strip())
    except FileNotFoundError:
        raise SystemExit(
            "ERROR: 'amrfinder' not found on PATH. Install AMRFinderPlus "
            "(conda install -c bioconda ncbi-amrfinderplus) and set up its "
            "database with `amrfinder -u`.")
    except subprocess.CalledProcessError as e:
        raise SystemExit(f"ERROR: `amrfinder --version` failed: {e.stderr or e}")


def run_amrfinder(fna, out_tsv, threads):
    """Run AMRFinderPlus on one nucleotide FASTA -> out_tsv. True on success."""
    cmd = ["amrfinder", "-n", str(fna), "-o", str(out_tsv), "--threads", str(threads)]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL,
                       stderr=subprocess.PIPE, text=True)
        return True
    except subprocess.CalledProcessError as e:
        out_tsv.unlink(missing_ok=True)  # drop partial output so a re-run retries
        last = e.stderr.strip().splitlines()[-1] if e.stderr else str(e)
        sys.stderr.write(f"WARNING: amrfinder failed on {fna.name}: {last}\n")
        return False


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assembly-dir", required=True, type=Path)
    ap.add_argument("--metadata", required=True, type=Path)
    ap.add_argument("--gene", action="append", required=True, metavar="NAME=PATTERN",
                    help="status column NAME from AMRFinderPlus gene-symbol regex "
                         "PATTERN (repeatable)")
    ap.add_argument("--amrfinder-out", type=Path, default=None,
                    help="dir for per-assembly AMRFinderPlus TSVs "
                         "(default: <assembly-dir>/.amrfinder)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--force", action="store_true",
                    help="re-run amrfinder even if a per-assembly TSV already exists")
    args = ap.parse_args(argv)

    gene_specs = parse_gene_specs(args.gene)
    print(f"AMRFinderPlus: {amrfinder_version()}")
    out_dir = args.amrfinder_out or (args.assembly_dir / ".amrfinder")
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Run AMRFinderPlus per assembly (resume-aware: skip existing TSVs).
    fnas = sorted(args.assembly_dir.glob("*.fna"))
    if not fnas:
        raise SystemExit(f"ERROR: no .fna files in {args.assembly_dir}")
    statuses = {}  # accession -> {column: positive|negative}
    for i, fna in enumerate(fnas, 1):
        acc = fna.stem
        tsv = out_dir / f"{acc}.tsv"
        if args.force or not tsv.exists():
            print(f"[{i}/{len(fnas)}] amrfinder {acc} ...", flush=True)
            if not run_amrfinder(fna, tsv, args.threads):
                continue  # leave unscored (blank)
        statuses[acc] = status_for(gene_symbols_from_tsv(tsv), gene_specs)

    # 2. Append status columns to metadata.csv, matched by accession (= .fna stem).
    if not args.metadata.exists():
        raise SystemExit(f"ERROR: metadata not found: {args.metadata}")
    with open(args.metadata, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])
    if "accession" not in fieldnames:
        raise SystemExit("ERROR: metadata.csv has no 'accession' column")
    for name, _ in gene_specs:
        if name not in fieldnames:
            fieldnames.append(name)

    n_scored = 0
    for row in rows:
        acc = row.get("accession", "")
        if acc in statuses:
            row.update(statuses[acc])
            n_scored += 1
        else:
            for name, _ in gene_specs:
                row.setdefault(name, "")  # unscored -> blank, never a silent negative

    with open(args.metadata, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    cols = ", ".join(n for n, _ in gene_specs)
    print(f"\nUpdated {n_scored}/{len(rows)} metadata rows with: {cols} -> {args.metadata}")
    for name, _ in gene_specs:
        pos = sum(1 for r in rows if r.get(name) == "positive")
        neg = sum(1 for r in rows if r.get(name) == "negative")
        blank = len(rows) - pos - neg
        print(f"  {name}: {pos} positive, {neg} negative, {blank} unscored/blank")
    return 0


if __name__ == "__main__":
    sys.exit(main())
