#!/usr/bin/env python3
"""host_category.py

OPTIONAL AssayVal companion helper — NOT part of the pipeline or the conda
environment, and not required by AssayVal.

Derive a `host_category` column (e.g. human / animal / environmental) for a
AssayVal metadata.csv from free-text host fields already present in the
metadata (by default `bs_host` and `bs_isolation_source`, captured by the
download helper). The column then becomes the grouping/scoring axis for a
host-associated marker assay, e.g. `target_group = host_category:human` for a
human fecal marker like HF183.

Host strings are unpredictable ("Homo sapiens", "cattle feces", "wastewater",
"river water", ...), so mapping is rule-based and inspectable:

1. See what's actually in your data first:
       python3 host_category.py --metadata metadata.csv --list-values
   This prints the distinct source-field values with counts, and exits.

2. Map with ordered, first-match-wins rules (CATEGORY=REGEX, case-insensitive,
   matched against the joined source fields). Provide your own with --rule,
   and/or start from the built-in set with --default-rules. Your --rule entries
   are checked BEFORE the built-in rules, so they override them:
       python3 host_category.py --metadata metadata.csv --default-rules
       python3 host_category.py --metadata metadata.csv \\
           --rule "human=wastewater|sewage" --default-rules

Unmatched rows get --default (blank by default), so nothing is silently
mis-categorized. Review the built-in rules (print them with --print-rules) —
the human/animal/environmental boundaries encode judgment calls (e.g. this set
classifies wastewater/sewage as ENVIRONMENTAL by sample matrix, not human; if
your source-tracking design treats wastewater as the human signal, add a
`--rule "human=wastewater|sewage"` to override).

Usage:
    python3 host_category.py --metadata metadata.csv --list-values
    python3 host_category.py --metadata metadata.csv --default-rules
"""
import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path

# Built-in starting rules (first match wins). Reviewed/overridable — see docstring.
DEFAULT_RULES = [
    ("human", r"homo sapiens|\bhuman\b|patient|clinical|nosocomial"),
    ("animal", r"bos taurus|cattle|\bcow\b|bovine|sus scrofa|swine|\bpig\b|porcine|"
               r"gallus|chicken|poultry|broiler|\bovis\b|sheep|goat|canis|\bdog\b|"
               r"canine|felis|\bcat\b|equus|horse|mus musculus|\bmouse\b|\brat\b|"
               r"rodent|avian|\bbird\b|livestock|\banimal\b"),
    ("environmental", r"wastewater|waste ?water|sewage|influent|effluent|\bwater\b|"
                      r"soil|sediment|marine|river|lake|ocean|\bsea\b|sludge|"
                      r"environment"),
]


def parse_rules(specs):
    """['category=regex', ...] -> [(category, compiled_regex), ...]."""
    out = []
    for s in specs:
        if "=" not in s:
            raise ValueError(f"--rule must be CATEGORY=REGEX (got '{s}')")
        cat, pattern = s.split("=", 1)
        cat = cat.strip()
        if not cat:
            raise ValueError(f"--rule has an empty category (got '{s}')")
        try:
            rx = re.compile(pattern, re.IGNORECASE)
        except re.error as e:
            raise ValueError(f"--rule '{cat}' has an invalid regex '{pattern}': {e}")
        out.append((cat, rx))
    return out


def compile_default_rules():
    return [(cat, re.compile(pat, re.IGNORECASE)) for cat, pat in DEFAULT_RULES]


def classify(text, rules, default=""):
    """Return the category of the first rule whose regex matches text, else default."""
    for cat, rx in rules:
        if rx.search(text):
            return cat
    return default


def _source_text(row, source_cols):
    return " ".join(str(row.get(c, "") or "") for c in source_cols)


def main(argv=None):
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", required=True, type=Path)
    ap.add_argument("--source-cols", default="bs_host,bs_isolation_source",
                    help="comma-separated metadata columns to read host text from "
                         "(default: bs_host,bs_isolation_source)")
    ap.add_argument("--column", default="host_category",
                    help="output column name (default: host_category)")
    ap.add_argument("--rule", action="append", default=[], metavar="CATEGORY=REGEX",
                    help="mapping rule, checked before the built-in rules (repeatable)")
    ap.add_argument("--default-rules", action="store_true",
                    help="append the built-in human/animal/environmental rules")
    ap.add_argument("--print-rules", action="store_true",
                    help="print the built-in rules and exit")
    ap.add_argument("--default", dest="default_cat", default="",
                    help="category for rows matching no rule (default: blank)")
    ap.add_argument("--list-values", action="store_true",
                    help="print distinct source-field values with counts, then exit")
    args = ap.parse_args(argv)

    if args.print_rules:
        for cat, pat in DEFAULT_RULES:
            print(f"{cat}={pat}")
        return 0

    source_cols = [c.strip() for c in args.source_cols.split(",") if c.strip()]
    if not args.metadata.exists():
        raise SystemExit(f"ERROR: metadata not found: {args.metadata}")
    with open(args.metadata, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
        fieldnames = list(reader.fieldnames or [])

    missing = [c for c in source_cols if c not in fieldnames]
    if missing:
        sys.stderr.write(f"WARNING: source column(s) not in metadata, ignored: {missing}\n")
    source_cols = [c for c in source_cols if c in fieldnames]
    if not source_cols:
        raise SystemExit(f"ERROR: none of the source columns exist in {args.metadata}")

    # --list-values: show what's actually there so rules can be written informedly.
    if args.list_values:
        counts = Counter(_source_text(r, source_cols).strip() for r in rows)
        print(f"# distinct values across {source_cols} ({len(counts)} unique):")
        for value, n in counts.most_common():
            print(f"{n:6d}  {value or '(empty)'}")
        return 0

    rules = parse_rules(args.rule)
    if args.default_rules:
        rules += compile_default_rules()
    if not rules:
        raise SystemExit(
            "ERROR: no rules given. Inspect your data with --list-values, then "
            "provide --rule CATEGORY=REGEX and/or --default-rules.")

    if args.column not in fieldnames:
        fieldnames.append(args.column)
    for row in rows:
        row[args.column] = classify(_source_text(row, source_cols), rules, args.default_cat)

    with open(args.metadata, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    tally = Counter(r[args.column] for r in rows)
    print(f"Wrote '{args.column}' for {len(rows)} rows -> {args.metadata}")
    for cat, n in tally.most_common():
        print(f"  {cat or '(unmatched/blank)'}: {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
