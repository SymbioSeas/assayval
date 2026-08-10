"""primeval command-line wrapper around the Snakemake workflow."""
import argparse
import re
import subprocess
import sys
import time
from pathlib import Path

from .rundir import resolve_run_dir

# <repo>/primeval/cli.py -> <repo>/workflow/Snakefile
WORKFLOW = Path(__file__).resolve().parent.parent / "workflow" / "Snakefile"

_PROGRESS_RE = re.compile(r'(\d+) of (\d+) steps \(([\d.]+)%\) done')
_ALERT_RE = re.compile(r'\b(error|exception|traceback|warning|failed)\b', re.IGNORECASE)


def filter_terminal_line(line: str, state: dict) -> list:
    """Return terminal lines for one raw Snakemake line (full line still logged)."""
    m = _PROGRESS_RE.search(line)
    if m:
        done, total, pct = m.group(1), m.group(2), float(m.group(3))
        milestone = int(pct // 10) * 10
        if milestone > state.get('last_milestone', -1):
            state['last_milestone'] = milestone
            return [f"  progress: {done}/{total} steps ({pct:.0f}%)"]
        return []
    if _ALERT_RE.search(line):
        return ["  " + line.rstrip()]
    return []


def build_parser():
    p = argparse.ArgumentParser(
        prog="primeval",
        description="Run the primeval in silico PCR pipeline on a set of assemblies.",
    )
    p.add_argument("--run-name", default="results",
                   help="Name for this run's results directory (default: results).")
    p.add_argument("--directory", default=".",
                   help="Analysis directory holding assemblies/, config.yaml, and outputs "
                        "(default: current directory).")
    p.add_argument("--configfile", default=None,
                   help="Pipeline config (default: <directory>/config.yaml).")
    p.add_argument("--force", action="store_true",
                   help="Reuse the existing dated run directory and resume unfinished work.")
    p.add_argument("--cores", type=int, default=8,
                   help="CPU cores for Snakemake (default: 8).")
    p.add_argument("snakemake_args", nargs=argparse.REMAINDER,
                   help="Arguments after -- are passed through to Snakemake.")
    return p


def build_snakemake_cmd(*, snakefile, workdir, configfile, results_dir_rel,
                        cores, quiet=True, extra=None):
    cmd = [
        "snakemake",
        "--snakefile", str(snakefile),
        "--directory", str(workdir),
        "--configfile", str(configfile),
        "--config", f"results_dir={results_dir_rel}",
        "--cores", str(cores),
    ]
    if quiet:
        cmd += ["--quiet", "rules"]
    cmd += list(extra or [])
    return cmd


def _load_config(configfile) -> dict:
    try:
        import yaml
        return yaml.safe_load(Path(configfile).read_text()) or {}
    except Exception:
        return {}


def _count_assemblies(configfile, workdir):
    try:
        cfg = _load_config(configfile)
        adir = Path(workdir) / cfg.get("assembly_dir", "assemblies")
        return len(list(adir.glob("*.fna")))
    except Exception:
        return None


DRY_RUN_FLAGS = {"-n", "--dry-run", "--dryrun"}


def _is_dry_run(extra):
    return any(a in DRY_RUN_FLAGS for a in extra)


def _tee_run(cmd, logf):
    lf = open(logf, "a") if logf else None
    state = {}
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in proc.stdout:
            if lf:
                lf.write(line)
            for out in filter_terminal_line(line, state):
                print(out, flush=True)
        proc.wait()
    finally:
        if lf:
            lf.close()
    return proc.returncode


def _print_header(run_name, rel_results, n_assemblies):
    n = "unknown" if n_assemblies is None else str(n_assemblies)
    print(f"primeval  |  run: {Path(rel_results).name}")
    print(f"  input   : assemblies/  ({n} assemblies)")
    print(f"  results : {rel_results}/")
    print()


def _print_footer(rel_results, elapsed):
    mins, secs = divmod(int(elapsed), 60)
    print(f"\nDone in {mins}m{secs:02d}s. Results in {rel_results}/")
    print(f"  {rel_results}/reports/assay_performance.csv       (per-assay sensitivity/specificity)")
    print(f"  {rel_results}/reports/detection_summary_long.csv  (per assay x group, all groupings)")
    print(f"  {rel_results}/reports/detection_by_assembly.csv   (per-assembly calls + metadata)")
    print(f"  {rel_results}/reports/figures/                    (one heatmap per grouping column)")


def main(argv=None, runner=None):
    args = build_parser().parse_args(argv)
    workdir = Path(args.directory).resolve()
    configfile = Path(args.configfile).resolve() if args.configfile else workdir / "config.yaml"
    if not configfile.exists():
        print(f"error: config file not found: {configfile}", file=sys.stderr)
        print(f"  copy the template:  cp <primeval-repo>/config/config.yaml "
              f"{workdir}/config.yaml", file=sys.stderr)
        return 2

    extra = list(args.snakemake_args)
    if extra and extra[0] == "--":
        extra = extra[1:]
    dry = _is_dry_run(extra)

    # The results root honors results_dir from the config (default "results");
    # the dated run directory is created inside it.
    results_root = workdir / _load_config(configfile).get("results_dir", "results")
    run_dir, mode = resolve_run_dir(results_root, args.run_name, args.force)
    try:
        rel_results = run_dir.relative_to(workdir).as_posix()
    except ValueError:  # absolute results_dir outside the analysis directory
        rel_results = run_dir.as_posix()
    # A dry run must not create the run directory or a run.log (it would burn the
    # run name and litter empty dated dirs). Only real runs materialize the dir.
    if dry:
        logf = None
    else:
        run_dir.mkdir(parents=True, exist_ok=True)
        logf = run_dir / "run.log"

    _print_header(args.run_name, rel_results, _count_assemblies(configfile, workdir))
    if not dry:
        print(f"Running primeval (full log: {rel_results}/run.log) ...")

    cmd = build_snakemake_cmd(
        snakefile=WORKFLOW, workdir=workdir, configfile=configfile,
        results_dir_rel=rel_results, cores=args.cores, extra=extra,
    )
    runner = runner or _tee_run
    t0 = time.perf_counter()
    rc = runner(cmd, logf)
    elapsed = time.perf_counter() - t0
    if rc != 0:
        print(f"\nerror: run failed (exit {rc}); see {logf}", file=sys.stderr)
    elif dry:
        print(f"\nDry run complete (no files written). Real run would write to {rel_results}/")
    else:
        _print_footer(rel_results, elapsed)
    return rc


if __name__ == "__main__":
    sys.exit(main())
