"""Resolve the results directory for an AssayVal run."""
from datetime import date as _date
from pathlib import Path


def resolve_run_dir(results_root, run_name, force, today=None, exists=None):
    """Return (run_dir, mode).

    mode is one of:
      "new"          -> results_root/<run_name>_<date> did not exist
      "resume"       -> it existed and force=True (reuse; Snakemake resumes)
      "incremented"  -> it existed and force=False (first free _<n>, n>=2)
    """
    results_root = Path(results_root)
    today = today or _date.today()
    exists = exists or (lambda p: Path(p).exists())
    datestr = today.strftime("%Y-%m-%d")
    base = results_root / f"{run_name}_{datestr}"
    if not exists(base):
        return base, "new"
    if force:
        return base, "resume"
    n = 2
    while exists(results_root / f"{run_name}_{datestr}_{n}"):
        n += 1
    return results_root / f"{run_name}_{datestr}_{n}", "incremented"
