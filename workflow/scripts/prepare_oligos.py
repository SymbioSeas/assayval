import re
import argparse
from pathlib import Path
from datetime import datetime

import table_io

# Handles slash-delimited (/ZEN/, /56-FAM/) and bracket-format ([AmMC6]) IDT tags.
# Phosphorothioate (*) notation is not stripped — notify if encountered in input.
_MOD_RE = re.compile(r'/[^/]+/|\[[^\]]+\]')


def strip_idt_modifications(seq: str) -> str:
    return _MOD_RE.sub('', seq)


def load_assay_table(csv_path: str) -> list[dict]:
    # Encoding/delimiter handling and validation live in table_io so every
    # reader of the assay table agrees on what it contains.
    return table_io.load_assay_table(csv_path)


def write_oligo_fasta(assays: list[dict], fasta_path: str, log_path: str) -> None:
    log_lines = [
        "# Oligo preparation log",
        f"# Generated: {datetime.now().isoformat()}",
    ]
    fasta_lines = []
    cleaned_oligos = []  # list of (oligo_id, cleaned_seq) for degenerate base scan
    for row in assays:
        name = row['assay']
        for role in ('fwd', 'rev', 'probe'):
            original = row[role]
            cleaned = strip_idt_modifications(original)
            oligo_id = f"{name}_{role}"
            if not cleaned:
                if role == 'probe':
                    log_lines.append(f"{oligo_id}: probe-free assay (SYBR/dsDNA dye) — no probe sequence written")
                else:
                    log_lines.append(f"{oligo_id}: skipped (empty sequence after stripping)")
                continue
            fasta_lines.append(f">{oligo_id}\n{cleaned}")
            cleaned_oligos.append((oligo_id, cleaned))
            if cleaned != original:
                log_lines.append(
                    f"{oligo_id}: stripped modification | original: {original} | cleaned: {cleaned}"
                )
            else:
                log_lines.append(f"{oligo_id}: no modification")
    Path(fasta_path).write_text('\n'.join(fasta_lines) + '\n')
    # Scan for IUPAC degenerate bases (beyond standard ACGTN)
    degenerate_found = False
    for oligo_id, seq in cleaned_oligos:
        degen = sorted(set(seq.upper()) & set('RYSWKMBDHV'))
        if degen:
            log_lines.append(
                f"IUPAC degenerate bases detected in {oligo_id}: {degen}"
            )
            degenerate_found = True
    if degenerate_found:
        log_lines.append("IUPAC-aware mismatch counting will be applied automatically.")
    else:
        log_lines.append("No degenerate bases detected in any oligo.")
    Path(log_path).write_text('\n'.join(log_lines) + '\n')


def main():
    parser = argparse.ArgumentParser(description="Prepare oligo FASTA from assay table")
    parser.add_argument('--assay-table', required=True)
    parser.add_argument('--fasta-out', required=True)
    parser.add_argument('--log-out', required=True)
    args = parser.parse_args()
    assays = load_assay_table(args.assay_table)
    write_oligo_fasta(assays, args.fasta_out, args.log_out)
    n_written = sum(1 for l in Path(args.fasta_out).read_text().splitlines() if l.startswith('>'))
    print(f"Wrote {n_written} oligos to {args.fasta_out}")


if __name__ == '__main__':
    main()
