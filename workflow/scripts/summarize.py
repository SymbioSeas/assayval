import table_io
import hashlib
import subprocess
import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from datetime import datetime


DETECTION_COLS = [
    'accession', 'assay', 'detection_call', 'n_amplicons', 'multi_amplicon_flag',
    'amplicon_sizes', 'contig_ids', 'fwd_mismatches', 'rev_mismatches',
    'probe_mismatches', 'probe_strand',
]

DETECTION_CALL_ORDER = ['Detected', 'Primer Only', 'Not Detected']

_HIGH_ANI_STATUSES = {'species_match', 'subspecies_match', 'derived_species_match'}


def _compute_ani_confidence(row) -> str:
    if row['ani_match_status'] in _HIGH_ANI_STATUSES and row['ani_taxonomy_check'] == 'OK':
        return 'High'
    if row['ani_match_status'] == 'genus_match' and row['ani_taxonomy_check'] == 'OK':
        return 'Genus'
    return 'Low'


def _effective_species(row) -> str:
    # genus_match assemblies fold into the single low-confidence bucket rather than
    # a separate per-species '(genus match)' group (which doubled the species axis of
    # every report). ani_confidence still records 'Genus' per-assembly, so the
    # detail survives in detection_by_assembly.csv and the manifest tier counts.
    if row['ani_confidence'] == 'High':
        return row['ani_best_match_organism']
    return 'Unclassified (low confidence ANI)'


def load_detection_results(amplicons_dir: str) -> pd.DataFrame:
    """Load per-assembly detection CSVs from a directory, skipping amplicons and empty files."""
    paths = [p for p in Path(amplicons_dir).glob('*.csv') if p.stat().st_size > 0]
    if not paths:
        return pd.DataFrame(columns=DETECTION_COLS)
    return pd.concat([pd.read_csv(p) for p in paths], ignore_index=True)


def join_metadata(det: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """Left-join detection results with ALL metadata columns on 'accession'.
    Only 'accession' is required in meta."""
    if 'accession' not in meta.columns:
        raise ValueError("metadata.csv must contain an 'accession' column")
    meta_cols = [c for c in meta.columns if c == 'accession' or c not in det.columns]
    return det.merge(meta[meta_cols], on='accession', how='left')


_ANI_COLS = {'ani_match_status', 'ani_taxonomy_check', 'ani_best_match_organism'}


def normalize_group_by(value) -> list:
    if not value:
        return []
    if isinstance(value, str):
        return [value]
    return list(value)


def _grouped_values(series: pd.Series) -> pd.Series:
    s = series.astype('object')
    return s.where(s.notna(), 'Ungrouped').replace('', 'Ungrouped')


def compute_grouping(joined: pd.DataFrame, group_by: list) -> tuple[pd.DataFrame, list]:
    """Return (joined, grouping_cols). Explicit columns, else ANI auto-grouping."""
    if group_by:
        missing = [c for c in group_by if c not in joined.columns]
        if missing:
            raise ValueError(
                f"group_by column(s) {missing} not found in metadata. "
                f"Available: {sorted(joined.columns)}")
        return joined, list(group_by)
    if _ANI_COLS.issubset(joined.columns):
        joined = joined.copy()
        joined['ani_best_match_organism'] = joined['ani_best_match_organism'].fillna('Unclassified')
        joined['ani_match_status'] = joined['ani_match_status'].fillna('no_match')
        joined['ani_taxonomy_check'] = joined['ani_taxonomy_check'].fillna('unknown')
        joined['ani_confidence'] = joined.apply(_compute_ani_confidence, axis=1)
        joined['group'] = joined.apply(_effective_species, axis=1)
        return joined, ['group']
    raise ValueError(
        "No group_by set and metadata lacks NCBI ANI columns. "
        "Set 'group_by' in config.yaml to the metadata column(s) to group by.")


GROUP_SUMMARY_COLS = [
    'grouping', 'group', 'assay', 'n_assemblies', 'n_detected',
    'n_primer_only', 'n_not_detected', 'pct_detected', 'pct_detected_or_primer',
    'n_multi_amplicon', 'max_amplicons',
]

DETECTION_SUMMARY_LONG_COLS = [
    'assay', 'grouping', 'group', 'n_assemblies', 'n_detected', 'n_primer_only',
    'n_not_detected', 'pct_detected', 'pct_detected_or_primer',
    'n_multi_amplicon', 'max_amplicons',
]


def build_group_summary(joined: pd.DataFrame, gcol: str) -> pd.DataFrame:
    """Compute per group-value × per-assay detection counts and percentages
    for one grouping column."""
    df = joined.copy()
    df['group'] = _grouped_values(df[gcol])
    rows = []
    for (grpname, assay), grp in df.groupby(['group', 'assay']):
        n = len(grp)
        n_det = int((grp['detection_call'] == 'Detected').sum())
        n_po = int((grp['detection_call'] == 'Primer Only').sum())
        n_nd = int((grp['detection_call'] == 'Not Detected').sum())
        namp = pd.to_numeric(grp['n_amplicons'], errors='coerce').fillna(0)
        rows.append({
            'grouping': gcol, 'group': grpname, 'assay': assay, 'n_assemblies': n,
            'n_detected': n_det, 'n_primer_only': n_po, 'n_not_detected': n_nd,
            'pct_detected': round(100 * n_det / n, 2) if n else 0.0,
            'pct_detected_or_primer': round(100 * (n_det + n_po) / n, 2) if n else 0.0,
            'n_multi_amplicon': int((namp > 1).sum()),
            'max_amplicons': int(namp.max()) if len(namp) else 0,
        })
    return pd.DataFrame(rows, columns=GROUP_SUMMARY_COLS)


def build_assay_summary_long(summary_all: pd.DataFrame) -> pd.DataFrame:
    """Per assay × grouping × group, ALL groups (including pct_detected=0),
    sorted by assay, grouping, then pct_detected descending."""
    if summary_all.empty:
        return pd.DataFrame(columns=DETECTION_SUMMARY_LONG_COLS)
    df = summary_all[DETECTION_SUMMARY_LONG_COLS].copy()
    return df.sort_values(['assay', 'grouping', 'pct_detected'],
                          ascending=[True, True, False]).reset_index(drop=True)


def build_group_matrix(summary: pd.DataFrame) -> pd.DataFrame:
    """Wide group × assay matrix of pct_detected with a leading n_assemblies column.
    Operates on a single-column summary (one grouping column)."""
    if summary.empty:
        return pd.DataFrame(columns=['group', 'n_assemblies'])
    n_assemblies = summary.groupby('group')['n_assemblies'].first()
    matrix = summary.pivot_table(index='group', columns='assay',
                                 values='pct_detected', aggfunc='first')
    matrix.insert(0, 'n_assemblies', n_assemblies)
    return matrix.reset_index()


def build_detection_matrix(joined: pd.DataFrame) -> pd.DataFrame:
    """Binary detection matrix: rows=accession, cols=assay, 1=Detected else 0."""
    dups = joined.duplicated(subset=['accession', 'assay'])
    if dups.any():
        raise ValueError(
            f"Duplicate accession×assay rows before pivot: "
            f"{joined[dups][['accession', 'assay']].values.tolist()}"
        )
    verbose = joined.pivot_table(
        index='accession', columns='assay', values='detection_call', aggfunc='first'
    )
    return verbose.map(lambda x: 1 if x == 'Detected' else 0)


def build_detection_by_assembly(binary: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    """One row per analyzed accession: all metadata columns + one 0/1 column per assay."""
    return meta.merge(binary.reset_index(), on='accession', how='right')


def write_group_heatmap(matrix: pd.DataFrame, figures_dir: str, gcol: str) -> None:
    """Group × assay heatmap colored by pct_detected (0-100)."""
    Path(figures_dir).mkdir(parents=True, exist_ok=True)
    if matrix.empty:
        return
    mat = matrix.set_index('group').drop(columns=['n_assemblies'], errors='ignore')
    if mat.empty:
        return
    annotate = mat.shape[0] * mat.shape[1] <= 300
    fig, ax = plt.subplots(figsize=(max(6, len(mat.columns)), max(4, len(mat) * 0.4)))
    sns.heatmap(mat, ax=ax, cmap='viridis', vmin=0, vmax=100,
                linewidths=0.5, linecolor='white', annot=annotate, fmt='.0f',
                cbar_kws={'label': '% detected'})
    ax.set_title(f'{gcol} detection (% of assemblies detected)')
    ax.set_xlabel('Assay')
    ax.set_ylabel(gcol)
    plt.tight_layout()
    for ext in ('pdf', 'png'):
        fig.savefig(f"{figures_dir}/{gcol}_detection_heatmap.{ext}", dpi=150)
    plt.close(fig)


_DETECTION_OUTPUT_COLS = [
    'detection_call', 'n_amplicons', 'multi_amplicon_flag', 'amplicon_sizes',
    'contig_ids', 'fwd_mismatches', 'rev_mismatches', 'probe_mismatches', 'probe_strand',
]


def dynamic_per_assay_columns(joined: pd.DataFrame) -> list:
    reserved = set(_DETECTION_OUTPUT_COLS) | {'accession', 'assay'}
    meta_cols = [c for c in joined.columns if c not in reserved]
    ordered = ['accession'] + meta_cols + \
        [c for c in _DETECTION_OUTPUT_COLS if c in joined.columns]
    return [c for c in ordered if c in joined.columns]


# Core columns of assay_performance.csv. Precision (PPV) is deliberately not
# reported: it is prevalence-dependent, so it moves with dataset composition
# rather than with assay quality. It stays recoverable as tp / (tp + fp).
ASSAY_PERF_COLS = [
    'assay', 'target_group', 'target_column', 'n_target', 'n_nontarget',
    'tp', 'fn', 'fp', 'tn', 'sensitivity', 'specificity',
    'n_detected_total',
]

# Assay-definition columns carried over from assay_table.csv into
# assay_performance.csv (each only when present), so the results file is
# self-contained and needs no cross-referencing back to the assay table.
ASSAY_CONTEXT_COLS = ['fwd', 'rev', 'probe', 'target_gene', *table_io.THRESHOLD_COLS]


def load_assay_targets(assay_table_path: str) -> dict:
    """Map assay name -> raw target_group string ('' when column absent or blank)."""
    rows = table_io.read_csv_rows(assay_table_path)
    return {r['assay']: (r.get('target_group') or '').strip() for r in rows}


def load_assay_context(assay_table_path: str) -> tuple[dict, list]:
    """Return ({assay: {col: value}}, present_cols) for the ASSAY_CONTEXT_COLS
    that exist in the assay table. Absent columns are simply not carried."""
    rows = table_io.read_csv_rows(assay_table_path)
    if not rows:
        return {}, []
    present = [c for c in ASSAY_CONTEXT_COLS if c in rows[0]]
    ctx = {r['assay']: {c: (r.get(c) or '') for c in present} for r in rows}
    return ctx, present


def _perf_columns(context_cols) -> list:
    """Column order: assay | carried assay-definition columns | target + metrics."""
    return ['assay'] + list(context_cols or []) + \
        [c for c in ASSAY_PERF_COLS if c != 'assay']


def parse_target(target: str, primary_col: str):
    """Return (column, value). 'col:val' -> (col, val); bare 'val' -> (primary_col, val);
    blank -> (None, None). Splits on the first colon only."""
    target = (target or '').strip()
    if not target:
        return None, None
    if ':' in target:
        col, val = target.split(':', 1)
        return col.strip(), val.strip()
    return primary_col, target


def _pct(num: int, den: int):
    return round(100 * num / den, 2) if den else None


def build_assay_performance(joined: pd.DataFrame, assay_targets: dict,
                            primary_col: str, assay_context: dict | None = None,
                            context_cols: list | None = None) -> pd.DataFrame:
    """Per-assay sensitivity/specificity against each assay's target.

    target_group is 'column:value' (or bare value -> primary_col). Positive
    signal = detection_call == 'Detected'. Assays without a target get a roster
    row with blank metrics. Missing target column -> ValueError; a value matching
    no assemblies -> warning. assay_context/context_cols carry assay-definition
    columns (fwd/rev/probe/target_gene) over from the assay table.
    """
    rows = []
    for assay, grp in joined.groupby('assay'):
        raw = (assay_targets.get(assay) or '').strip()
        ctx = (assay_context or {}).get(assay, {})
        detected = grp['detection_call'] == 'Detected'
        n_detected_total = int(detected.sum())
        col, val = parse_target(raw, primary_col)
        if col is None:
            rows.append({
                'assay': assay, **ctx, 'target_group': '', 'target_column': '',
                'n_target': None, 'n_nontarget': None, 'tp': None, 'fn': None,
                'fp': None, 'tn': None, 'sensitivity': None, 'specificity': None,
                'n_detected_total': n_detected_total,
            })
            continue
        if col not in grp.columns:
            raise ValueError(
                f"assay '{assay}' target references column '{col}' not in metadata")
        is_target = grp[col].astype('object') == val
        n_target = int(is_target.sum())
        if n_target == 0:
            print(f"WARNING: assay '{assay}' target '{raw}' matched 0 assemblies "
                  f"(check spelling / that target genomes are in your dataset).")
        n_nontarget = int((~is_target).sum())
        tp = int((is_target & detected).sum())
        fn = n_target - tp
        fp = int((~is_target & detected).sum())
        tn = n_nontarget - fp
        rows.append({
            'assay': assay, **ctx, 'target_group': raw, 'target_column': col,
            'n_target': n_target, 'n_nontarget': n_nontarget,
            'tp': tp, 'fn': fn, 'fp': fp, 'tn': tn,
            'sensitivity': _pct(tp, n_target),
            'specificity': _pct(tn, n_nontarget),
            'n_detected_total': n_detected_total,
        })
    df = pd.DataFrame(rows, columns=_perf_columns(context_cols))
    _count_cols = ['n_target', 'n_nontarget', 'tp', 'fn', 'fp', 'tn', 'n_detected_total']
    df[_count_cols] = df[_count_cols].astype('Int64')
    return df


def _md5(path: str) -> str:
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(65536), b''):
            h.update(chunk)
    return h.hexdigest()


def _tool_version(cmd: list) -> str:
    try:
        return subprocess.check_output(cmd, stderr=subprocess.STDOUT,
                                       text=True).strip().split('\n')[0]
    except Exception:
        return 'unknown'


def _per_assay_settings(assay_table: str) -> list:
    """Manifest lines for assays whose own threshold columns override the
    config, and for every oligo carrying LNA bases (1-based positions)."""
    from oligo import parse_oligo, OligoError
    lines = []
    for r in table_io.read_csv_rows(assay_table):
        name = r.get('assay', '')
        overrides = [f"{c}={r[c]}" for c in table_io.THRESHOLD_COLS if (r.get(c) or '').strip()]
        if overrides:
            lines.append(f"{name}: {', '.join(overrides)}")
        for col in ('fwd', 'rev', 'probe'):
            # Validation is enforced upstream (prepare_oligos/detect); the
            # manifest only reports, so an unparsable cell is skipped here.
            try:
                lna = parse_oligo(r.get(col, '')).lna
            except OligoError:
                continue
            if lna:
                lines.append(f"{name}: {col} LNA at {', '.join(str(i + 1) for i in lna)}")
    return lines


def write_run_manifest(manifest_path: str, params: dict, assay_table: str, counts: dict | None = None) -> None:
    """Write run manifest with parameters, tool versions, and checksums."""
    lines = [
        f"# Run manifest — generated {datetime.now().isoformat()}",
        "",
        "## Parameters",
    ]
    for k, v in params.items():
        lines.append(f"  {k}: {v}")
    lines += [
        "",
        "## Tool versions",
        f"  BLAST: {_tool_version(['blastn', '-version'])}",
        f"  Python: {_tool_version(['python', '--version'])}",
        "",
        "## Checksums",
        f"  assay_table.csv  MD5: {_md5(assay_table)}",
    ]
    per_assay = _per_assay_settings(assay_table)
    if per_assay:
        lines += ["", "## Per-assay settings"] + [f"  {l}" for l in per_assay]
    if counts:
        lines += ["", "## Input accounting"]
        for k, v in counts.items():
            lines.append(f"  {k}: {v}")
    Path(manifest_path).write_text('\n'.join(lines) + '\n')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--amplicons-dir', required=True)
    p.add_argument('--metadata', required=True)
    p.add_argument('--assay-table', required=True)
    p.add_argument('--reports-dir', required=True)
    p.add_argument('--group-by', nargs='*', default=[])
    p.add_argument('--max-primer-mismatches', type=int, default=2)
    p.add_argument('--prime3-exact-nt', type=int, default=1)
    p.add_argument('--max-probe-mismatches', type=int, default=1)
    p.add_argument('--max-amplicon-size', type=int, default=500)
    p.add_argument('--store-amplicon-sequences',
                   type=lambda x: x.lower() == 'true', default=True)
    p.add_argument('--keep-blast', type=lambda x: x.lower() == 'true', default=False)
    p.add_argument('--keep-logs', type=lambda x: x.lower() == 'true', default=False)
    p.add_argument('--amplicon-fasta', type=lambda x: x.lower() == 'true', default=True)
    p.add_argument('--amplicon-flank-bp', type=int, default=50)
    p.add_argument('--lna-mismatch', default='exact')
    p.add_argument('--blast-params', default='',
                   help="Comma-separated key=value BLAST search parameters, recorded "
                        "in run_manifest.txt. Empty on a re-scored run, which inherits "
                        "its BLAST output from the run named by --rescored-from.")
    p.add_argument('--rescored-from', default='',
                   help="Path of the run whose cached BLAST output this run re-scored "
                        "(empty on a normal run); recorded in run_manifest.txt.")
    args = p.parse_args()

    reports = Path(args.reports_dir)
    figures = reports / 'figures'
    per_assay = reports / 'per_assay'
    per_assay.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    det = load_detection_results(args.amplicons_dir)
    meta = table_io.read_metadata(args.metadata)
    joined = join_metadata(det, meta)
    joined, grouping_cols = compute_grouping(joined, list(args.group_by))
    primary_col = grouping_cols[0]

    # Per-assay tables (dynamic columns)
    cols = dynamic_per_assay_columns(joined)
    for assay, grp in joined.groupby('assay'):
        sorted_grp = grp.assign(
            _call_order=pd.Categorical(
                grp['detection_call'], categories=DETECTION_CALL_ORDER, ordered=True
            )
        ).sort_values('_call_order')[cols]
        sorted_grp.to_csv(per_assay / f"{assay}_results.csv", index=False)

    # One matrix + heatmap per grouping column; combined long table
    summaries = []
    for gcol in grouping_cols:
        summ = build_group_summary(joined, gcol)
        summaries.append(summ)
        matrix = build_group_matrix(summ)
        matrix.to_csv(reports / f"{gcol}_detection_matrix.csv", index=False)
        write_group_heatmap(matrix, str(figures), gcol)
    summary_all = pd.concat(summaries, ignore_index=True) if summaries \
        else pd.DataFrame(columns=GROUP_SUMMARY_COLS)

    detection_long = build_assay_summary_long(summary_all)
    detection_long.to_csv(reports / 'detection_summary_long.csv', index=False)
    with pd.ExcelWriter(reports / 'assay_summary.xlsx', engine='openpyxl') as writer:
        for assay, grp in detection_long.groupby('assay'):
            grp.to_excel(writer, sheet_name=str(assay)[:31], index=False)

    # Per-assay sensitivity / specificity
    assay_targets = load_assay_targets(args.assay_table)
    assay_context, context_cols = load_assay_context(args.assay_table)
    build_assay_performance(joined, assay_targets, primary_col,
                            assay_context=assay_context,
                            context_cols=context_cols).to_csv(
        reports / 'assay_performance.csv', index=False)

    # Per-assembly detection joined with input metadata
    binary = build_detection_matrix(joined)
    build_detection_by_assembly(binary, meta).to_csv(
        reports / 'detection_by_assembly.csv', index=False)

    # Run manifest
    params = {
        'max_primer_mismatches': args.max_primer_mismatches,
        'prime3_exact_nt': args.prime3_exact_nt,
        'max_probe_mismatches': args.max_probe_mismatches,
        'max_amplicon_size': args.max_amplicon_size,
        'store_amplicon_sequences': args.store_amplicon_sequences,
        'keep_blast': args.keep_blast,
        'keep_logs': args.keep_logs,
        'amplicon_fasta': args.amplicon_fasta,
        'amplicon_flank_bp': args.amplicon_flank_bp,
        'lna_mismatch': args.lna_mismatch,
    }
    # BLAST search provenance. Present only on a run that performed the search;
    # a re-scored run records rescored_from instead, pointing at the manifest
    # that does hold them.
    for kv in filter(None, (s.strip() for s in args.blast_params.split(','))):
        key, _, value = kv.partition('=')
        params[f'blast_{key}'] = value
    # Provenance for re-scored runs: which run's cached BLAST output this one
    # was derived from. Omitted entirely on a normal run so manifests from
    # ordinary runs are unchanged.
    if args.rescored_from:
        params['rescored_from'] = args.rescored_from
    counts = {
        'assemblies_analyzed': int(joined['accession'].nunique()),
        'metadata_rows': int(len(meta)),
        'grouping_columns': ", ".join(grouping_cols),
        'n_assays': int(joined['assay'].nunique()),
    }
    # In ANI-auto grouping mode, record how many assemblies were species-confident
    # vs genus-only vs low-confidence (genus/low share the 'Unclassified' group).
    if 'ani_confidence' in joined.columns:
        conf = joined.drop_duplicates('accession')['ani_confidence'].value_counts()
        counts['ani_high_confidence'] = int(conf.get('High', 0))
        counts['ani_genus_only'] = int(conf.get('Genus', 0))
        counts['ani_low_confidence'] = int(conf.get('Low', 0))
    write_run_manifest(str(reports / 'run_manifest.txt'), params, args.assay_table, counts=counts)

    print(f"Reports written to {args.reports_dir}")


if __name__ == '__main__':
    main()
