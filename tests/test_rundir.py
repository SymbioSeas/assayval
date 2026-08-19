from datetime import date
from pathlib import Path
from assayval.rundir import resolve_run_dir

D = date(2026, 7, 7)


def _exists_set(paths):
    s = {str(p) for p in paths}
    return lambda p: str(p) in s


def test_new_when_absent():
    root = Path("/wd/results")
    run_dir, mode = resolve_run_dir(root, "Vpop", force=False, today=D, exists=_exists_set([]))
    assert run_dir == root / "Vpop_2026-07-07"
    assert mode == "new"


def test_resume_when_exists_and_force():
    root = Path("/wd/results")
    base = root / "Vpop_2026-07-07"
    run_dir, mode = resolve_run_dir(root, "Vpop", force=True, today=D, exists=_exists_set([base]))
    assert run_dir == base
    assert mode == "resume"


def test_increment_when_exists_no_force():
    root = Path("/wd/results")
    base = root / "Vpop_2026-07-07"
    run_dir, mode = resolve_run_dir(root, "Vpop", force=False, today=D, exists=_exists_set([base]))
    assert run_dir == root / "Vpop_2026-07-07_2"
    assert mode == "incremented"


def test_increment_skips_existing_suffixes():
    root = Path("/wd/results")
    existing = [root / "Vpop_2026-07-07", root / "Vpop_2026-07-07_2"]
    run_dir, mode = resolve_run_dir(root, "Vpop", force=False, today=D, exists=_exists_set(existing))
    assert run_dir == root / "Vpop_2026-07-07_3"


def test_default_run_name_results():
    root = Path("/wd/results")
    run_dir, _ = resolve_run_dir(root, "results", force=False, today=D, exists=_exists_set([]))
    assert run_dir == root / "results_2026-07-07"
