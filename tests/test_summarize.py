import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))

import pytest
import pandas as pd
import numpy as np
from summarize import (
    load_detection_results,
    join_metadata,
    normalize_group_by,
    compute_grouping,
    build_group_summary,
    build_group_matrix,
    build_assay_summary_long,
    build_detection_matrix,
    build_detection_by_assembly,
    write_group_heatmap,
    write_run_manifest,
    load_assay_targets,
    parse_target,
    build_assay_performance,
    dynamic_per_assay_columns,
    ASSAY_PERF_COLS,
)


# --- Fixtures ---

def make_detection_df(accessions, assays, call='Detected', n_amplicons=1):
    """Build a minimal detection DataFrame (like output of run_detection)."""
    rows = []
    for acc in accessions:
        for assay in assays:
            rows.append({
                'accession': acc, 'assay': assay, 'detection_call': call,
                'n_amplicons': n_amplicons, 'multi_amplicon_flag': n_amplicons > 1,
                'amplicon_sizes': '150', 'contig_ids': 'c1',
                'fwd_mismatches': 0, 'rev_mismatches': 0,
                'probe_mismatches': 0, 'probe_strand': '+',
            })
    return pd.DataFrame(rows)


def make_metadata_df(accessions, species='Vibrio harveyi',
                     ani_match_status='species_match', ani_taxonomy_check='OK'):
    rows = []
    for acc in accessions:
        rows.append({
            'accession': acc,
            'ani_best_match_organism': species,
            'ani_match_status': ani_match_status,
            'ani_taxonomy_check': ani_taxonomy_check,
            'infraspecific_strain': 'strain_X',
            'bs_strain': 'bs_X',
            'bs_host': 'shrimp',
            'bs_isolation_source': 'water',
            'assembly_level': 'Complete Genome',
        })
    return pd.DataFrame(rows)


# --- load_detection_results ---

def test_load_detection_results(tmp_path):
    """Load multiple per-assembly CSVs into one DataFrame."""
    assays = ['VhPath', 'Valg']
    for i, acc in enumerate(['GCF_001', 'GCF_002']):
        df = make_detection_df([acc], assays)
        df.to_csv(tmp_path / f"{acc}.csv", index=False)
    combined = load_detection_results(str(tmp_path))
    assert len(combined) == 4  # 2 assemblies × 2 assays
    assert set(combined['accession']) == {'GCF_001', 'GCF_002'}


def test_load_detection_results_empty_dir(tmp_path):
    """Empty directory returns empty DataFrame."""
    df = load_detection_results(str(tmp_path))
    assert df.empty


def test_load_detection_results_skips_empty_files(tmp_path):
    """Zero-byte CSV files are skipped without error."""
    (tmp_path / "empty.csv").write_text("")
    df = make_detection_df(['GCF_001'], ['VhPath'])
    df.to_csv(tmp_path / "GCF_001.csv", index=False)
    combined = load_detection_results(str(tmp_path))
    assert len(combined) == 1
    assert combined.iloc[0]['accession'] == 'GCF_001'


# --- join_metadata ---

def test_join_metadata_basic():
    det = make_detection_df(['GCF_001'], ['VhPath'])
    meta = make_metadata_df(['GCF_001'])
    result = join_metadata(det, meta)
    assert 'ani_best_match_organism' in result.columns
    assert result.iloc[0]['ani_best_match_organism'] == 'Vibrio harveyi'
    assert result.iloc[0]['bs_host'] == 'shrimp'


def test_join_metadata_unclassified():
    """Assemblies with no metadata match get NaN metadata from the plain join;
    ANI auto-grouping (compute_grouping) resolves this to 'Unclassified' / Low."""
    det = make_detection_df(['GCF_999'], ['VhPath'])
    meta = make_metadata_df(['GCF_001'])  # GCF_999 not in metadata
    joined, _ = compute_grouping(join_metadata(det, meta), [])
    assert joined.iloc[0]['ani_best_match_organism'] == 'Unclassified'
    assert joined.iloc[0]['ani_confidence'] == 'Low'


def test_join_metadata_minimal_only_accession():
    det = make_detection_df(['GCF_001'], ['A'])
    meta = pd.DataFrame([{'accession': 'GCF_001', 'phenotype': 'protective'}])
    joined = join_metadata(det, meta)
    assert joined.iloc[0]['phenotype'] == 'protective'


# --- normalize_group_by / compute_grouping ---

def test_normalize_group_by():
    assert normalize_group_by(None) == []
    assert normalize_group_by("") == []
    assert normalize_group_by("species") == ["species"]
    assert normalize_group_by(["species", "phenotype"]) == ["species", "phenotype"]


def test_compute_grouping_explicit_multi_column():
    det = make_detection_df(['GCF_001', 'GCF_002'], ['A'])
    meta = pd.DataFrame([
        {'accession': 'GCF_001', 'species': 'Vibrio mediterranei', 'phenotype': 'protective'},
        {'accession': 'GCF_002', 'species': 'Vibrio harveyi', 'phenotype': ''},
    ])
    joined, cols = compute_grouping(join_metadata(det, meta), ['species', 'phenotype'])
    assert cols == ['species', 'phenotype']


def test_compute_grouping_missing_column_raises():
    det = make_detection_df(['GCF_001'], ['A'])
    meta = pd.DataFrame([{'accession': 'GCF_001', 'phenotype': 'x'}])
    with pytest.raises(ValueError, match="group_by"):
        compute_grouping(join_metadata(det, meta), ['not_a_column'])


def test_compute_grouping_ani_auto_mode():
    det = make_detection_df(['GCF_001'], ['A'])
    meta = make_metadata_df(['GCF_001'], species='Vibrio harveyi')
    joined, cols = compute_grouping(join_metadata(det, meta), [])
    assert cols == ['group']
    assert joined.iloc[0]['group'] == 'Vibrio harveyi'
    assert joined.iloc[0]['ani_confidence'] == 'High'


def test_compute_grouping_no_group_by_no_ani_raises():
    det = make_detection_df(['GCF_001'], ['A'])
    meta = pd.DataFrame([{'accession': 'GCF_001', 'phenotype': 'x'}])
    with pytest.raises(ValueError, match="group_by"):
        compute_grouping(join_metadata(det, meta), [])


# --- build_assay_summary_long / build_group_matrix ---

def test_build_assay_summary_long_includes_misses_and_sorts():
    det = pd.concat([
        make_detection_df(['GCF_001', 'GCF_002'], ['VhPath'], call='Detected'),
        make_detection_df(['GCF_003'], ['VhPath'], call='Not Detected', n_amplicons=0),
    ])
    meta = pd.concat([
        make_metadata_df(['GCF_001', 'GCF_002'], species='Vibrio harveyi'),
        make_metadata_df(['GCF_003'], species='Vibrio sp.'),
    ])
    joined, _ = compute_grouping(join_metadata(det, meta), [])
    summary = build_group_summary(joined, 'group')
    long = build_assay_summary_long(summary)
    # Miss (pct_detected=0) is INCLUDED
    assert (long['pct_detected'] == 0).any()
    assert list(long.columns) == [
        'assay', 'grouping', 'group', 'n_assemblies', 'n_detected', 'n_primer_only',
        'n_not_detected', 'pct_detected', 'pct_detected_or_primer',
        'n_multi_amplicon', 'max_amplicons',
    ]
    # Sorted by assay, then pct_detected desc (hit before miss within VhPath)
    vh = long[long['assay'] == 'VhPath']['pct_detected'].tolist()
    assert vh == sorted(vh, reverse=True)


def test_build_group_matrix_is_wide():
    det = pd.concat([
        make_detection_df(['GCF_001'], ['VhPath'], call='Detected'),
        make_detection_df(['GCF_001'], ['Valg'], call='Not Detected', n_amplicons=0),
    ])
    meta = make_metadata_df(['GCF_001'])
    joined, _ = compute_grouping(join_metadata(det, meta), [])
    matrix = build_group_matrix(build_group_summary(joined, 'group'))
    assert 'group' in matrix.columns
    assert 'n_assemblies' in matrix.columns
    assert 'VhPath' in matrix.columns and 'Valg' in matrix.columns
    r = matrix[matrix['group'] == 'Vibrio harveyi'].iloc[0]
    assert r['VhPath'] == 100.0 and r['Valg'] == 0.0


# --- build_group_summary ---

def test_build_group_summary_counts():
    det = pd.concat([
        make_detection_df(['GCF_001', 'GCF_002'], ['VhPath'], call='Detected'),
        make_detection_df(['GCF_003'], ['VhPath'], call='Not Detected'),
        make_detection_df(['GCF_004'], ['VhPath'], call='Primer Only'),
    ])
    meta = make_metadata_df(['GCF_001', 'GCF_002', 'GCF_003', 'GCF_004'])
    joined, _ = compute_grouping(join_metadata(det, meta), [])
    summary = build_group_summary(joined, 'group')
    row = summary[(summary['group'] == 'Vibrio harveyi') &
                  (summary['assay'] == 'VhPath')].iloc[0]
    assert row['n_assemblies'] == 4
    assert row['n_detected'] == 2
    assert row['n_not_detected'] == 1
    assert row['n_primer_only'] == 1
    assert abs(row['pct_detected'] - 50.0) < 0.01
    assert abs(row['pct_detected_or_primer'] - 75.0) < 0.01
    assert row['n_multi_amplicon'] == 0
    assert row['max_amplicons'] == 1


def test_build_group_summary_empty_cell_is_ungrouped():
    det = make_detection_df(['GCF_001', 'GCF_002'], ['A'])
    meta = pd.DataFrame([
        {'accession': 'GCF_001', 'phenotype': 'protective'},
        {'accession': 'GCF_002', 'phenotype': ''},
    ])
    joined, _ = compute_grouping(join_metadata(det, meta), ['phenotype'])
    summ = build_group_summary(joined, 'phenotype')
    assert set(summ['grouping']) == {'phenotype'}
    assert 'Ungrouped' in set(summ['group'])


def test_accession_absent_from_metadata_is_ungrouped():
    det = make_detection_df(['GCF_999'], ['A'])
    meta = pd.DataFrame([{'accession': 'GCF_001', 'phenotype': 'x'}])
    joined, _ = compute_grouping(join_metadata(det, meta), ['phenotype'])
    summ = build_group_summary(joined, 'phenotype')
    assert summ.iloc[0]['group'] == 'Ungrouped'


# --- per-assay sort order ---

def test_per_assay_sort_order(tmp_path):
    """Per-assay CSVs are sorted Detected → Primer Only → Not Detected."""
    from summarize import DETECTION_CALL_ORDER

    # Build joined df with all three call types in scrambled order
    det = pd.concat([
        make_detection_df(['GCF_003'], ['VhPath'], call='Not Detected'),
        make_detection_df(['GCF_002'], ['VhPath'], call='Primer Only'),
        make_detection_df(['GCF_001'], ['VhPath'], call='Detected'),
    ])
    meta = make_metadata_df(['GCF_001', 'GCF_002', 'GCF_003'])
    joined = join_metadata(det, meta)

    per_assay_dir = tmp_path / "per_assay"
    per_assay_dir.mkdir()
    cols = dynamic_per_assay_columns(joined)
    for assay, grp in joined.groupby('assay'):
        sorted_grp = grp.assign(
            _call_order=pd.Categorical(
                grp['detection_call'], categories=DETECTION_CALL_ORDER, ordered=True
            )
        ).sort_values('_call_order')[cols]
        sorted_grp.to_csv(per_assay_dir / f"{assay}_results.csv", index=False)

    result = pd.read_csv(per_assay_dir / "VhPath_results.csv")
    calls = result['detection_call'].tolist()
    assert calls == ['Detected', 'Primer Only', 'Not Detected']


# --- build_detection_matrix ---

def test_build_detection_matrix_binary():
    det = pd.concat([
        make_detection_df(['GCF_001'], ['VhPath'], call='Detected'),
        make_detection_df(['GCF_001'], ['Valg'], call='Not Detected', n_amplicons=0),
    ])
    joined = join_metadata(det, make_metadata_df(['GCF_001']))
    binary = build_detection_matrix(joined)
    assert binary.loc['GCF_001', 'VhPath'] == 1
    assert binary.loc['GCF_001', 'Valg'] == 0


def test_build_detection_matrix_duplicate_raises():
    det = pd.concat([
        make_detection_df(['GCF_001'], ['VhPath'], call='Detected'),
        make_detection_df(['GCF_001'], ['VhPath'], call='Not Detected'),
    ])
    joined = join_metadata(det, make_metadata_df(['GCF_001']))
    with pytest.raises(ValueError, match="Duplicate"):
        build_detection_matrix(joined)


def test_build_detection_by_assembly_joins_metadata():
    det = pd.concat([
        make_detection_df(['GCF_001'], ['VhPath'], call='Detected'),
        make_detection_df(['GCF_001'], ['Valg'], call='Not Detected', n_amplicons=0),
    ])
    meta = make_metadata_df(['GCF_001'])
    joined = join_metadata(det, meta)
    binary = build_detection_matrix(joined)
    out = build_detection_by_assembly(binary, meta)
    assert len(out) == 1
    row = out.iloc[0]
    assert row['accession'] == 'GCF_001'
    assert row['bs_host'] == 'shrimp'         # metadata carried through
    assert row['VhPath'] == 1 and row['Valg'] == 0


def test_write_group_heatmap_creates_files(tmp_path):
    det = pd.concat([
        make_detection_df(['GCF_001'], ['VhPath'], call='Detected'),
        make_detection_df(['GCF_002'], ['VhPath'], call='Not Detected', n_amplicons=0),
        make_detection_df(['GCF_001', 'GCF_002'], ['Valg'], call='Detected'),
    ])
    meta = make_metadata_df(['GCF_001', 'GCF_002'])
    joined, _ = compute_grouping(join_metadata(det, meta), [])
    matrix = build_group_matrix(build_group_summary(joined, 'group'))
    figures_dir = str(tmp_path / "figures")
    write_group_heatmap(matrix, figures_dir, 'group')
    assert (tmp_path / "figures" / "group_detection_heatmap.pdf").exists()
    assert (tmp_path / "figures" / "group_detection_heatmap.png").exists()


# --- write_run_manifest ---

def test_write_run_manifest(tmp_path):
    params = {
        'max_primer_mismatches': 2,
        'prime3_exact_nt': 3,
        'max_probe_mismatches': 1,
        'max_amplicon_size': 500,
        'store_amplicon_sequences': True,
    }
    assay_table_path = tmp_path / "assay_table.csv"
    assay_table_path.write_text("assay,probe,fwd,rev\n")
    write_run_manifest(str(tmp_path / "manifest.txt"), params, str(assay_table_path))
    content = (tmp_path / "manifest.txt").read_text()
    assert 'max_primer_mismatches: 2' in content
    assert 'prime3_exact_nt: 3' in content
    assert 'assay_table.csv' in content
    assert 'md5:' in content.lower() or 'MD5' in content
    assert 'blast' in content.lower() or 'BLAST' in content


def test_manifest_includes_input_accounting(tmp_path):
    at = tmp_path / "assay_table.csv"; at.write_text("assay,probe,fwd,rev\n")
    counts = {'assemblies_analyzed': 10, 'metadata_rows': 12,
              'grouping_columns': 'species, phenotype', 'n_assays': 5}
    write_run_manifest(str(tmp_path / "m.txt"), {'max_amplicon_size': 500},
                       str(at), counts=counts)
    txt = (tmp_path / "m.txt").read_text()
    assert "Input accounting" in txt
    assert "assemblies_analyzed: 10" in txt
    assert "grouping_columns: species, phenotype" in txt


# --- ANI confidence tiers ---

def test_ani_confidence_high():
    """species_match + OK → High."""
    det = make_detection_df(['GCF_001'], ['VhPath'])
    meta = make_metadata_df(['GCF_001'], ani_match_status='species_match', ani_taxonomy_check='OK')
    joined, _ = compute_grouping(join_metadata(det, meta), [])
    assert joined.iloc[0]['ani_confidence'] == 'High'


def test_ani_confidence_genus():
    """genus_match + OK → Genus."""
    det = make_detection_df(['GCF_001'], ['VhPath'])
    meta = make_metadata_df(['GCF_001'], ani_match_status='genus_match', ani_taxonomy_check='OK')
    joined, _ = compute_grouping(join_metadata(det, meta), [])
    assert joined.iloc[0]['ani_confidence'] == 'Genus'


def test_ani_confidence_low_mismatch():
    """mismatch + Inconclusive → Low."""
    det = make_detection_df(['GCF_001'], ['VhPath'])
    meta = make_metadata_df(['GCF_001'], ani_match_status='mismatch', ani_taxonomy_check='Inconclusive')
    joined, _ = compute_grouping(join_metadata(det, meta), [])
    assert joined.iloc[0]['ani_confidence'] == 'Low'


def test_group_summary_folds_genus_into_unclassified():
    """Genus-only assemblies fold into the single 'Unclassified' bucket, not a separate
    '(genus match)' group; ani_confidence still records 'Genus' per-assembly."""
    det = pd.concat([
        make_detection_df(['GCF_001'], ['VhPath'], call='Detected'),
        make_detection_df(['GCF_002'], ['VhPath'], call='Detected'),
    ])
    meta = pd.concat([
        make_metadata_df(['GCF_001'], ani_match_status='species_match', ani_taxonomy_check='OK'),
        make_metadata_df(['GCF_002'], ani_match_status='genus_match', ani_taxonomy_check='OK'),
    ])
    joined, _ = compute_grouping(join_metadata(det, meta), [])
    summary = build_group_summary(joined, 'group')
    groups = set(summary['group'])
    assert 'Vibrio harveyi' in groups
    assert 'Vibrio harveyi (genus match)' not in groups
    assert 'Unclassified (low confidence ANI)' in groups
    # per-assembly ani_confidence still distinguishes the genus-only assembly
    assert set(joined['ani_confidence']) == {'High', 'Genus'}
    high = summary[summary['group'] == 'Vibrio harveyi'].iloc[0]
    uncl = summary[summary['group'] == 'Unclassified (low confidence ANI)'].iloc[0]
    assert high['n_assemblies'] == 1
    assert uncl['n_assemblies'] == 1


# --- assay_performance.csv (column:value targets) ---

def test_load_assay_targets_reads_column(tmp_path):
    p = tmp_path / "assays.csv"
    p.write_text("assay,probe,fwd,rev,target_group\nA,,F,R,phenotype:protective\nB,,F,R,\n")
    assert load_assay_targets(str(p)) == {"A": "phenotype:protective", "B": ""}


def test_load_assay_targets_missing_column(tmp_path):
    p = tmp_path / "assays.csv"
    p.write_text("assay,probe,fwd,rev\nA,,F,R\n")
    assert load_assay_targets(str(p)) == {"A": ""}


def test_parse_target_forms():
    assert parse_target("phenotype:protective", "species") == ("phenotype", "protective")
    assert parse_target("Vibrio mediterranei", "species") == ("species", "Vibrio mediterranei")
    assert parse_target("", "species") == (None, None)
    # value may contain a colon; only the first splits
    assert parse_target("note:a:b", "species") == ("note", "a:b")


def test_build_assay_performance_column_value_target():
    # species-wide assay targets species:Vmed (all phenotypes); phenotype col also present
    rows = []
    data = [
        ("a1", "Vibrio mediterranei", "protective", "Detected"),
        ("a2", "Vibrio mediterranei", "pathogenic", "Detected"),
        ("a3", "Vibrio harveyi", "", "Detected"),      # off-target FP
        ("a4", "Vibrio harveyi", "", "Not Detected"),  # TN
    ]
    for acc, sp, ph, call in data:
        rows.append({"accession": acc, "assay": "Vmed", "species": sp,
                     "phenotype": ph, "detection_call": call})
    joined = pd.DataFrame(rows)
    perf = build_assay_performance(joined, {"Vmed": "species:Vibrio mediterranei"}, "species")
    r = perf.iloc[0]
    assert list(perf.columns) == ASSAY_PERF_COLS
    assert r["target_column"] == "species"
    assert r["n_target"] == 2 and r["n_nontarget"] == 2
    assert r["tp"] == 2 and r["fn"] == 0 and r["fp"] == 1 and r["tn"] == 1
    assert r["sensitivity"] == 100.0 and r["specificity"] == 50.0


def test_build_assay_performance_nested_phenotype_target():
    # clade assay targets phenotype:protective; intermediate/pathogenic Vmed are non-target
    rows = []
    data = [
        ("a1", "Vibrio mediterranei", "protective", "Detected"),   # TP
        ("a2", "Vibrio mediterranei", "protective", "Not Detected"),  # FN
        ("a3", "Vibrio mediterranei", "pathogenic", "Detected"),   # FP (wrong clade)
        ("a4", "Vibrio harveyi", "", "Not Detected"),              # TN
    ]
    for acc, sp, ph, call in data:
        rows.append({"accession": acc, "assay": "VmProt", "species": sp,
                     "phenotype": ph, "detection_call": call})
    joined = pd.DataFrame(rows)
    perf = build_assay_performance(joined, {"VmProt": "phenotype:protective"}, "species")
    r = perf.iloc[0]
    assert r["target_column"] == "phenotype"
    assert r["n_target"] == 2 and r["tp"] == 1 and r["fn"] == 1
    assert r["fp"] == 1 and r["tn"] == 1
    assert r["sensitivity"] == 50.0 and r["specificity"] == 50.0


def test_build_assay_performance_primer_only_is_negative():
    joined = pd.DataFrame([
        {"accession": "a1", "assay": "A", "species": "t", "detection_call": "Primer Only"},
        {"accession": "a2", "assay": "A", "species": "t", "detection_call": "Detected"},
    ])
    perf = build_assay_performance(joined, {"A": "t"}, "species")
    r = perf.iloc[0]
    assert r["tp"] == 1 and r["fn"] == 1 and r["sensitivity"] == 50.0


def test_build_assay_performance_no_target_blank_metrics():
    joined = pd.DataFrame([
        {"accession": "a1", "assay": "Ctrl", "species": "x", "detection_call": "Detected"},
    ])
    perf = build_assay_performance(joined, {"Ctrl": ""}, "species")
    r = perf.iloc[0]
    assert r["target_group"] == "" and r["target_column"] == "" and pd.isna(r["sensitivity"])
    assert r["n_detected_total"] == 1


def test_build_assay_performance_int_counts_with_blank_target_row():
    # Mixing a targeted assay with a blank-target (roster) assay used to upcast
    # the count columns to float (NaN from the blank row), rendering counts like
    # "1.0" in the CSV. Counts should stay pandas nullable Int64 so the targeted
    # row's counts are true integers and the blank row's counts are <NA>.
    joined = pd.DataFrame([
        {"accession": "a1", "assay": "A", "species": "t", "detection_call": "Detected"},
        {"accession": "a2", "assay": "A", "species": "x", "detection_call": "Not Detected"},
        {"accession": "a1", "assay": "Ctrl", "species": "t", "detection_call": "Detected"},
    ])
    perf = build_assay_performance(joined, {"A": "t", "Ctrl": ""}, "species")
    count_cols = ["n_target", "n_nontarget", "tp", "fn", "fp", "tn", "n_detected_total"]
    for col in count_cols:
        assert str(perf[col].dtype) == "Int64"

    targeted = perf[perf["assay"] == "A"].iloc[0]
    assert targeted["n_target"] == 1 and targeted["n_nontarget"] == 1
    assert targeted["tp"] == 1 and targeted["fn"] == 0
    assert targeted["fp"] == 0 and targeted["tn"] == 1

    blank = perf[perf["assay"] == "Ctrl"].iloc[0]
    assert pd.isna(blank["n_target"]) and pd.isna(blank["tp"])

    csv_text = perf.to_csv(index=False)
    assert ",1.0," not in csv_text
    assert not csv_text.rstrip("\n").endswith("1.0")


def test_build_assay_performance_missing_column_raises():
    joined = pd.DataFrame([
        {"accession": "a1", "assay": "A", "species": "t", "detection_call": "Detected"},
    ])
    with pytest.raises(ValueError, match="nope"):
        build_assay_performance(joined, {"A": "nope:x"}, "species")


def test_load_assay_context_present_and_absent_columns(tmp_path):
    """Only ASSAY_CONTEXT_COLS actually present in the assay table are carried."""
    from summarize import load_assay_context
    p = tmp_path / "assays.csv"
    p.write_text(
        "assay,probe,fwd,rev,target_gene,notes\n"
        "A,PPP,FFF,RRR,gyrB,some note\n")
    ctx, present = load_assay_context(str(p))
    # order follows ASSAY_CONTEXT_COLS; 'notes' is not carried
    assert present == ["fwd", "rev", "probe", "target_gene"]
    assert ctx["A"] == {"fwd": "FFF", "rev": "RRR", "probe": "PPP", "target_gene": "gyrB"}

    p2 = tmp_path / "assays2.csv"
    p2.write_text("assay,probe,fwd,rev\nA,,FFF,RRR\n")   # no target_gene
    ctx2, present2 = load_assay_context(str(p2))
    assert present2 == ["fwd", "rev", "probe"]
    assert ctx2["A"] == {"fwd": "FFF", "rev": "RRR", "probe": ""}


def test_build_assay_performance_carries_assay_context():
    """Assay-definition columns land between 'assay' and the target/metric columns."""
    joined = pd.DataFrame([
        {"accession": "a1", "assay": "A", "species": "t", "detection_call": "Detected"},
        {"accession": "a2", "assay": "Ctrl", "species": "t", "detection_call": "Detected"},
    ])
    ctx = {
        "A": {"fwd": "FFF", "rev": "RRR", "probe": "PPP", "target_gene": "gyrB"},
        "Ctrl": {"fwd": "F2", "rev": "R2", "probe": "", "target_gene": ""},
    }
    cols = ["fwd", "rev", "probe", "target_gene"]
    perf = build_assay_performance(
        joined, {"A": "species:t", "Ctrl": ""}, "species",
        assay_context=ctx, context_cols=cols)
    assert list(perf.columns) == [
        "assay", "fwd", "rev", "probe", "target_gene", "target_group", "target_column",
        "n_target", "n_nontarget", "tp", "fn", "fp", "tn",
        "sensitivity", "specificity", "n_detected_total",
    ]
    a = perf[perf["assay"] == "A"].iloc[0]
    assert a["fwd"] == "FFF" and a["probe"] == "PPP" and a["target_gene"] == "gyrB"
    assert a["sensitivity"] == 100.0
    # context is carried for blank-target control rows too
    c = perf[perf["assay"] == "Ctrl"].iloc[0]
    assert c["fwd"] == "F2" and pd.isna(c["sensitivity"])


def test_build_assay_performance_no_precision_column():
    """Precision (PPV) is deliberately not reported; tp/fp remain so it stays derivable."""
    joined = pd.DataFrame([
        {"accession": "a1", "assay": "A", "species": "t", "detection_call": "Detected"},
    ])
    perf = build_assay_performance(joined, {"A": "species:t"}, "species")
    assert "precision" not in perf.columns
    assert {"tp", "fp"} <= set(perf.columns)


def test_build_assay_performance_zero_match_target_warns(capsys):
    """A target that matches no assemblies emits a warning and NA metrics."""
    joined = pd.DataFrame([
        {"accession": "a1", "assay": "A", "species": "Vibrio harveyi", "detection_call": "Detected"},
    ])
    perf = build_assay_performance(joined, {"A": "species:Vibrio mediterranei"}, "species")
    r = perf.iloc[0]
    assert r["n_target"] == 0
    assert pd.isna(r["sensitivity"])
    assert "matched 0 assemblies" in capsys.readouterr().out


# --- main() end-to-end integration ---

def test_summarize_main_end_to_end_multi_group(tmp_path):
    import subprocess, sys
    from pathlib import Path as P
    amps = tmp_path / "amplicons"; amps.mkdir()
    make_detection_df(['GCF_001'], ['Vmed']).to_csv(amps / "GCF_001.csv", index=False)
    make_detection_df(['GCF_002'], ['Vmed'], call='Not Detected', n_amplicons=0).to_csv(
        amps / "GCF_002.csv", index=False)
    (tmp_path / "metadata.csv").write_text(
        "accession,species,phenotype\n"
        "GCF_001,Vibrio mediterranei,protective\n"
        "GCF_002,Vibrio harveyi,\n")
    (tmp_path / "assays.csv").write_text(
        "assay,probe,fwd,rev,target_group\nVmed,,F,R,species:Vibrio mediterranei\n")
    reports = tmp_path / "reports"
    script = P(__file__).parent.parent / "workflow" / "scripts" / "summarize.py"
    subprocess.run([sys.executable, str(script),
        "--amplicons-dir", str(amps), "--metadata", str(tmp_path / "metadata.csv"),
        "--assay-table", str(tmp_path / "assays.csv"), "--reports-dir", str(reports),
        "--group-by", "species", "phenotype"], check=True)
    for f in ["species_detection_matrix.csv", "phenotype_detection_matrix.csv",
              "detection_summary_long.csv", "assay_performance.csv",
              "detection_by_assembly.csv", "run_manifest.txt",
              "figures/species_detection_heatmap.pdf",
              "figures/phenotype_detection_heatmap.pdf"]:
        assert (reports / f).exists(), f
    perf = pd.read_csv(reports / "assay_performance.csv").iloc[0]
    assert perf["target_column"] == "species"
    assert perf["sensitivity"] == 100.0 and perf["specificity"] == 100.0
    long = pd.read_csv(reports / "detection_summary_long.csv")
    assert set(long["grouping"]) == {"species", "phenotype"}


def test_summarize_main_ani_auto_manifest_tier_counts(tmp_path):
    """ANI-auto mode (no --group-by): genus folds into Unclassified and the manifest
    records high/genus/low tier counts."""
    import subprocess, sys
    from pathlib import Path as P
    amps = tmp_path / "amplicons"; amps.mkdir()
    for acc in ['GCF_001', 'GCF_002', 'GCF_003']:
        make_detection_df([acc], ['Vmed']).to_csv(amps / f"{acc}.csv", index=False)
    (tmp_path / "metadata.csv").write_text(
        "accession,ani_best_match_organism,ani_match_status,ani_taxonomy_check\n"
        "GCF_001,Vibrio mediterranei,species_match,OK\n"
        "GCF_002,Vibrio mediterranei,genus_match,OK\n"
        "GCF_003,Vibrio sp.,mismatch,Inconclusive\n")
    (tmp_path / "assays.csv").write_text("assay,probe,fwd,rev\nVmed,,F,R\n")
    reports = tmp_path / "reports"
    script = P(__file__).parent.parent / "workflow" / "scripts" / "summarize.py"
    subprocess.run([sys.executable, str(script),
        "--amplicons-dir", str(amps), "--metadata", str(tmp_path / "metadata.csv"),
        "--assay-table", str(tmp_path / "assays.csv"), "--reports-dir", str(reports)],
        check=True)
    # ANI-auto → single 'group' matrix; genus folds into Unclassified (no separate genus group)
    groups = set(pd.read_csv(reports / "group_detection_matrix.csv")["group"])
    assert "Vibrio mediterranei" in groups
    assert not any("genus match" in g for g in groups)
    assert "Unclassified (low confidence ANI)" in groups
    manifest = (reports / "run_manifest.txt").read_text()
    assert "ani_high_confidence: 1" in manifest
    assert "ani_genus_only: 1" in manifest
    assert "ani_low_confidence: 1" in manifest
