"""Write one amplicon FASTA per assay from the per-assembly amplicon records.

Each record is one valid amplicon (Detected or Primer Only) emitted by
detect.py, oriented fwd primer -> rev primer, with lowercase flanking sequence
around the uppercase primer-to-primer amplicon. The files are meant for
downstream design work, e.g. importing into Geneious to place a probe on a
probe-free assay or to inspect primer-site variation in off-target taxa.

Run once per run by the export_amplicon_fasta Snakemake rule
(workflow/rules/export.smk).
"""
import argparse
import re
from pathlib import Path

import pandas as pd

import table_io
from detect import AMPLICON_RECORD_COLS

FASTA_LINE_WIDTH = 70
UNKNOWN_ORGANISM = 'unknown'


def sanitize(text: str) -> str:
    """Make text safe for a FASTA ID / file name: anything outside
    [A-Za-z0-9.-] becomes '_' (runs collapsed, ends trimmed)."""
    return re.sub(r'[^A-Za-z0-9.-]+', '_', str(text)).strip('_')


def organism_lookup(meta: pd.DataFrame, group_by: list) -> dict:
    """accession -> organism label for headers.

    organism_name (written by download-assemblies) is preferred; a bring-your-own
    metadata.csv without it falls back to the first group_by column, then
    'unknown'. Empty cells are 'unknown' too.
    """
    if 'organism_name' in meta.columns:
        col = 'organism_name'
    elif group_by and group_by[0] in meta.columns:
        col = group_by[0]
    else:
        col = None
    out = {}
    for _, row in meta.iterrows():
        value = row[col] if col else None
        out[row['accession']] = (str(value).strip()
                                 if value is not None and pd.notna(value) and str(value).strip()
                                 else UNKNOWN_ORGANISM)
    return out


def _fmt(value) -> str:
    if value is None or (not isinstance(value, str) and pd.isna(value)):
        return 'NA'
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def format_header(rec: dict, organism: str) -> tuple[str, str]:
    """Return (fasta_id, header_line) for one amplicon record.

    The ID — the part Geneious and most tools use as the sequence name — is
    accession-first so it is unique and searchable, then the copy number
    (amp{k}of{n}), then the organism so alignment rows show taxa.
    """
    fasta_id = (f"{rec['accession']}_amp{_fmt(rec['amplicon_index'])}"
                f"of{_fmt(rec['n_amplicons'])}_{sanitize(organism)}")
    fields = [
        f'organism="{organism}"',
        f"assay={rec['assay']}",
        f"assembly_call={str(rec['assembly_call']).replace(' ', '_')}",
        f"probe={rec['probe_status']}",
        f"contig={rec['contig_id']}:{_fmt(rec['amplicon_start'])}-"
        f"{_fmt(rec['amplicon_end'])}({rec['strand']})",
        f"amplicon_bp={_fmt(rec['amplicon_size'])}",
        f"flank_bp={_fmt(rec['flank_up_bp'])},{_fmt(rec['flank_down_bp'])}",
        f"fwd_mm={_fmt(rec['fwd_mismatches'])}",
        f"rev_mm={_fmt(rec['rev_mismatches'])}",
        f"probe_mm={_fmt(rec['probe_mismatches'])}",
    ]
    return fasta_id, f">{fasta_id} " + " ".join(fields)


def _wrap(seq: str) -> str:
    return "\n".join(seq[i:i + FASTA_LINE_WIDTH]
                     for i in range(0, len(seq), FASTA_LINE_WIDTH))


def export_amplicon_fasta(records: pd.DataFrame, assays: list, organisms: dict,
                          out_dir) -> dict:
    """Write {assay}_amplicons.fasta for every assay (empty when it has no
    amplicons) plus amplicon_index.csv; return {assay: n_records}."""
    out_dir = Path(out_dir)
    files = {a: f"{sanitize(a)}_amplicons.fasta" for a in assays}
    seen = {}
    for assay, fname in files.items():
        if fname in seen:
            raise ValueError(f"assays '{seen[fname]}' and '{assay}' map to the same "
                             f"file name {fname}; rename one in assay_table.csv")
        seen[fname] = assay
    out_dir.mkdir(parents=True, exist_ok=True)

    records = records.sort_values(['assay', 'accession', 'amplicon_index'], kind='stable')
    index_rows, counts = [], {}
    for assay in assays:
        grp = records[records['assay'] == assay]
        lines = []
        for rec in grp.to_dict('records'):
            organism = organisms.get(rec['accession'], UNKNOWN_ORGANISM)
            fasta_id, header = format_header(rec, organism)
            lines += [header, _wrap(rec['sequence'])]
            row = {'fasta_id': fasta_id, 'organism': organism, 'fasta_file': files[assay]}
            row.update({k: v for k, v in rec.items() if k != 'sequence'})
            index_rows.append(row)
        (out_dir / files[assay]).write_text("\n".join(lines) + ("\n" if lines else ""))
        counts[assay] = len(grp)

    index_cols = ['fasta_id', 'organism', 'fasta_file'] + \
        [c for c in records.columns if c != 'sequence']
    pd.DataFrame(index_rows, columns=index_cols).to_csv(
        out_dir / 'amplicon_index.csv', index=False)
    return counts


def load_records(records_dir) -> pd.DataFrame:
    frames = [pd.read_csv(p, sep='\t') for p in sorted(Path(records_dir).glob('*.tsv'))
              if p.stat().st_size > 0]
    frames = [f for f in frames if not f.empty]
    if not frames:
        return pd.DataFrame(columns=AMPLICON_RECORD_COLS)
    return pd.concat(frames, ignore_index=True)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--records-dir', required=True,
                   help="Directory of per-assembly amplicon record TSVs from detect.py.")
    p.add_argument('--metadata', required=True)
    p.add_argument('--assay-table', required=True)
    p.add_argument('--group-by', nargs='*', default=[])
    p.add_argument('--out-dir', required=True)
    args = p.parse_args()

    records = load_records(args.records_dir)
    assays = [r['assay'] for r in table_io.load_assay_table(args.assay_table)]
    organisms = organism_lookup(table_io.read_metadata(args.metadata), args.group_by)
    counts = export_amplicon_fasta(records, assays, organisms, args.out_dir)
    total = sum(counts.values())
    print(f"Wrote {total} amplicons across {len(counts)} assay FASTA files to {args.out_dir}")


if __name__ == '__main__':
    main()
