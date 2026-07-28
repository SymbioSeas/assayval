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

    Note: order-independence assumes unique 'accession' values across records;
    duplicate/missing accessions are ties in the sort and fall back to input order.
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
