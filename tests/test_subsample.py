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
