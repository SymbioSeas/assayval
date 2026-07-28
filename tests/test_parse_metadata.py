import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts" / "download"))

from parse_metadata import merge_rows


def test_merge_unions_by_accession_new_wins():
    existing = [{"accession": "A", "organism_name": "old"},
                {"accession": "B", "organism_name": "b"}]
    new = [{"accession": "A", "organism_name": "new"},
           {"accession": "C", "organism_name": "c"}]
    out = merge_rows(existing, new)
    by_acc = {r["accession"]: r["organism_name"] for r in out}
    assert by_acc == {"A": "new", "B": "b", "C": "c"}   # A updated, B kept, C added
    assert [r["accession"] for r in out] == ["A", "B", "C"]  # existing order, new appended


def test_merge_into_existing_csv_end_to_end(tmp_path):
    import subprocess, sys as _sys, json, csv as _csv
    out_csv = tmp_path / "metadata.csv"
    # seed an existing metadata.csv (run 1)
    j1 = tmp_path / "r1.jsonl"
    j1.write_text(json.dumps({"accession": "GCF_1", "organism": {"organism_name": "Taxon A"}}) + "\n")
    script = Path(__file__).parent.parent / "scripts" / "download" / "parse_metadata.py"
    subprocess.run([_sys.executable, str(script), "--jsonl", str(j1), "--out", str(out_csv), "--merge"], check=True)
    # run 2 into the same CSV with --merge
    j2 = tmp_path / "r2.jsonl"
    j2.write_text(json.dumps({"accession": "GCF_2", "organism": {"organism_name": "Taxon B"}}) + "\n")
    subprocess.run([_sys.executable, str(script), "--jsonl", str(j2), "--out", str(out_csv), "--merge"], check=True)
    accs = {row["accession"] for row in _csv.DictReader(open(out_csv))}
    assert accs == {"GCF_1", "GCF_2"}   # both runs present, not just the last
