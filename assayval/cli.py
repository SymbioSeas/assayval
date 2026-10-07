"""AssayVal command-line wrapper around the Snakemake workflow."""
import argparse
import re
import subprocess
import sys
import time
from pathlib import Path

from .rundir import resolve_run_dir

# <repo>/assayval/cli.py -> <repo>/workflow/Snakefile
WORKFLOW = Path(__file__).resolve().parent.parent / "workflow" / "Snakefile"
REPO_ROOT = Path(__file__).resolve().parent.parent

# Starter files `assay-val init` drops into an analysis directory, as
# (destination name, source in the repo). config.yaml is copied straight from
# the repo default so a new project always starts from the current defaults.
INIT_TEMPLATES = [
    ("config.yaml", REPO_ROOT / "config" / "config.yaml"),
    ("assay_table.csv", REPO_ROOT / "config" / "assay_table_template.csv"),
]


def version_string() -> str:
    """Package version, plus the git commit when running from a checkout.

    An analysis is only reproducible if the exact code can be recovered, and a
    release version alone does not identify a working tree. When AssayVal is
    installed editable from a clone (the normal case for a pipeline run), the
    short commit and a -dirty marker are appended so a run log pins the code.
    """
    try:
        from importlib.metadata import version as _pkg_version
        ver = _pkg_version("assayval")
    except Exception:
        ver = "unknown"
    try:
        sha = subprocess.check_output(
            ["git", "-C", str(REPO_ROOT), "rev-parse", "--short", "HEAD"],
            text=True, stderr=subprocess.DEVNULL).strip()
        dirty = subprocess.run(
            ["git", "-C", str(REPO_ROOT), "diff", "--quiet"],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode != 0
        ver += f" ({sha}{'-dirty' if dirty else ''})"
    except Exception:
        pass
    return ver

_PROGRESS_RE = re.compile(r'(\d+) of (\d+) steps \(([\d.]+)%\) done')
# "<path>.py:<line>: SomeWarning: <message>" as printed by Python's warnings module.
_PYWARN_RE = re.compile(r'^\S+\.py:\d+: \w*Warning: (.*)$')
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
    m = _PYWARN_RE.match(line.strip())
    if m:
        # Python warnings from the workflow scripts (e.g. an assay table read as
        # Mac Roman) repeat once per job; show each distinct one once.
        msg = f"  warning: {m.group(1)}"
        seen = state.setdefault('seen_warnings', set())
        if msg in seen:
            return []
        seen.add(msg)
        return [msg]
    if _ALERT_RE.search(line):
        return ["  " + line.rstrip()]
    return []


def build_parser():
    p = argparse.ArgumentParser(
        prog="assay-val",
        description="Run the AssayVal in silico PCR pipeline on a set of assemblies.",
        epilog="To start a new analysis, run `assay-val init` in an empty directory to "
               "drop in a starter config.yaml and assay_table.csv "
               "(see `assay-val init --help`).",
    )
    p.add_argument("--version", action="version",
                   version=f"assay-val {version_string()}",
                   help="Print the version (and git commit, when run from a clone) and exit.")
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
    p.add_argument("--rescore-from", default=None, metavar="RUN_DIR",
                   help="Re-score a previous run's retained BLAST output instead of running "
                        "BLAST again. Point at that run's directory (or its blast/ directory). "
                        "The source run must have used keep_blast: true, and its assemblies "
                        "must still be present. Combine with --set to vary detection "
                        "thresholds cheaply, e.g. for a threshold sensitivity analysis.")
    p.add_argument("--set", dest="overrides", action="append", default=[],
                   metavar="KEY=VALUE",
                   help="Override one config value for this run; repeatable. Mainly for "
                        "detection thresholds, e.g. --set max_primer_mismatches=1 "
                        "--set prime3_exact_nt=2.")
    p.add_argument("snakemake_args", nargs=argparse.REMAINDER,
                   help="Arguments after -- are passed through to Snakemake.")
    return p


def build_init_parser():
    p = argparse.ArgumentParser(
        prog="assay-val init",
        description="Copy the starter config.yaml and assay_table.csv into an analysis "
                    "directory. Edit both for your project, then run `assay-val` from "
                    "that directory.",
    )
    p.add_argument("--directory", default=".",
                   help="Analysis directory to set up; created if missing "
                        "(default: current directory).")
    p.add_argument("--force", action="store_true",
                   help="Overwrite existing config.yaml / assay_table.csv.")
    return p


def init_main(argv=None):
    """Drop the starter templates into an analysis directory.

    Existing files are never overwritten without --force: they may already hold
    a project's edited thresholds or assays. Missing files are still written,
    and the exit code is 1 if anything was skipped so scripts notice.
    """
    args = build_init_parser().parse_args(argv)
    workdir = Path(args.directory).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    print(f"assay-val init  |  {workdir}")
    skipped = []
    for name, src in INIT_TEMPLATES:
        dest = workdir / name
        if dest.exists() and not args.force:
            skipped.append(name)
            print(f"  skipped : {name}  (already exists)")
            continue
        verb = "replaced" if dest.exists() else "wrote"
        dest.write_bytes(src.read_bytes())
        print(f"  {verb:<8}: {name}")

    if skipped:
        print(f"\nerror: left {len(skipped)} existing file(s) untouched; "
              f"re-run with --force to overwrite", file=sys.stderr)
        return 1

    print("\nNext steps:")
    print("  1. Edit assay_table.csv: replace the example rows with your assays.")
    print("  2. Edit config.yaml: check assembly_dir, metadata, thresholds, group_by.")
    print("  3. Get assemblies, e.g.:  download-assemblies -t \"<taxon>\" -o assemblies/")
    print("  4. Run from this directory:  assay-val --run-name <name>")
    return 0


def parse_overrides(overrides):
    """Validate --set KEY=VALUE tokens and return them unchanged.

    Values are passed through to Snakemake's --config verbatim, which parses them
    as YAML, so ints/bools/strings all behave as they do in config.yaml.
    """
    for token in overrides:
        key, sep, _ = token.partition("=")
        if not sep or not key.strip():
            raise ValueError(f"--set expects KEY=VALUE, got: {token!r}")
    return list(overrides)


def override_value(overrides, key):
    """Return the effective --set value for `key`, or None.

    Later flags win, matching how Snakemake resolves repeated --config keys, so
    the CLI reports the value the workflow will actually use.
    """
    value = None
    for token in overrides:
        k, sep, v = token.partition("=")
        if sep and k.strip() == key:
            value = v
    return value


def resolve_rescore_blast_dir(workdir, rescore_from):
    """Resolve --rescore-from to (blast_dir_for_config, n_cached_tsv).

    Accepts either a previous run directory or that run's blast/ directory
    directly. Raises ValueError with an actionable message when the cache is
    absent or empty — the usual cause is a source run left at the default
    keep_blast: false, which deletes each BLAST TSV as it is consumed.
    """
    workdir = Path(workdir)
    src = Path(rescore_from)
    if not src.is_absolute():
        src = workdir / src
    src = src.resolve()
    if not src.exists():
        raise ValueError(f"--rescore-from path does not exist: {src}")
    blast_dir = src if src.name == "blast" else src / "blast"
    if not blast_dir.is_dir():
        raise ValueError(
            f"no blast/ directory in {src}\n"
            "  --rescore-from needs a run that retained its raw BLAST output.\n"
            "  Re-run the source analysis with keep_blast: true in config.yaml."
        )
    n = len(list(blast_dir.glob("*.tsv.gz")))
    if n == 0:
        raise ValueError(
            f"no cached BLAST output (*.tsv.gz) in {blast_dir}\n"
            "  Re-run the source analysis with keep_blast: true in config.yaml."
        )
    try:
        as_config = blast_dir.relative_to(workdir).as_posix()
    except ValueError:  # cache lives outside the analysis directory
        as_config = blast_dir.as_posix()
    return as_config, n


def build_snakemake_cmd(*, snakefile, workdir, configfile, results_dir_rel,
                        cores, quiet=True, extra=None, config_overrides=None):
    # results_dir and any overrides share a single --config: Snakemake takes
    # KEY=VALUE with nargs='+', and a second --config flag would replace the
    # first rather than extend it.
    cmd = [
        "snakemake",
        "--snakefile", str(snakefile),
        "--directory", str(workdir),
        "--configfile", str(configfile),
        "--config", f"results_dir={results_dir_rel}",
        *list(config_overrides or []),
        "--cores", str(cores),
    ]
    if quiet:
        cmd += ["--quiet", "rules"]
    cmd += list(extra or [])
    return cmd


def _load_config(configfile) -> dict:
    try:
        import yaml
        # Explicit UTF-8 (not the locale default); undecodable bytes can only be
        # in comments here, so replace rather than fail the header/run setup.
        text = Path(configfile).read_text(encoding="utf-8", errors="replace")
        return yaml.safe_load(text) or {}
    except Exception:
        return {}


def _count_assemblies(configfile, workdir, assembly_dir=None):
    """Count input assemblies, honouring a --set assembly_dir override.

    Without this the header reports whatever the config file says while the run
    uses the override — which is exactly the kind of quiet mismatch a run log is
    supposed to prevent.
    """
    try:
        cfg = _load_config(configfile)
        adir = Path(workdir) / (assembly_dir or cfg.get("assembly_dir", "assemblies"))
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


def _print_header(run_name, rel_results, n_assemblies, rescore=None,
                  overrides=None, assembly_dir=None):
    print(f"assay-val  |  run: {Path(rel_results).name}")
    if rescore:
        blast_dir, n_cached = rescore
        print(f"  mode    : re-score (BLAST skipped)")
        print(f"  input   : {blast_dir}  ({n_cached} cached BLAST files)")
    else:
        n = "unknown" if n_assemblies is None else str(n_assemblies)
        print(f"  input   : {assembly_dir or 'assemblies'}/  ({n} assemblies)")
    if overrides:
        print(f"  config  : {' '.join(overrides)}")
    print(f"  results : {rel_results}/")
    print()


def _print_footer(rel_results, elapsed, run_dir=None):
    mins, secs = divmod(int(elapsed), 60)
    print(f"\nDone in {mins}m{secs:02d}s. Results in {rel_results}/")
    print(f"  {rel_results}/reports/assay_performance.csv       (per-assay sensitivity/specificity)")
    print(f"  {rel_results}/reports/detection_summary_long.csv  (per assay x group, all groupings)")
    print(f"  {rel_results}/reports/detection_by_assembly.csv   (per-assembly calls + metadata)")
    print(f"  {rel_results}/reports/figures/                    (one heatmap per grouping column)")
    if run_dir is not None and (Path(run_dir) / "amplicon_fasta").is_dir():
        print(f"  {rel_results}/amplicon_fasta/                     (one amplicon FASTA per assay)")


def main(argv=None, runner=None):
    argv = sys.argv[1:] if argv is None else list(argv)
    # `init` is dispatched ahead of the run parser so the long-standing flat run
    # interface (assay-val --run-name ...) is unchanged.
    if argv and argv[0] == "init":
        return init_main(argv[1:])

    args = build_parser().parse_args(argv)
    workdir = Path(args.directory).resolve()
    configfile = Path(args.configfile).resolve() if args.configfile else workdir / "config.yaml"
    if not configfile.exists():
        print(f"error: config file not found: {configfile}", file=sys.stderr)
        print(f"  create starter files here with:  assay-val init --directory {workdir}",
              file=sys.stderr)
        return 2

    try:
        user_overrides = parse_overrides(args.overrides)
    except ValueError as e:
        print(f"error: {e}", file=sys.stderr)
        return 2
    overrides = list(user_overrides)

    rescore = None
    if args.rescore_from:
        try:
            blast_dir, n_cached = resolve_rescore_blast_dir(workdir, args.rescore_from)
        except ValueError as e:
            print(f"error: {e}", file=sys.stderr)
            return 2
        rescore = (blast_dir, n_cached)
        # blast_dir switches the workflow into re-score mode (see workflow/Snakefile).
        overrides = overrides + [f"blast_dir={blast_dir}"]

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

    assembly_dir = override_value(user_overrides, "assembly_dir")
    _print_header(args.run_name, rel_results,
                  _count_assemblies(configfile, workdir, assembly_dir),
                  rescore=rescore, overrides=user_overrides,
                  assembly_dir=assembly_dir)
    if not dry:
        print(f"Running AssayVal (full log: {rel_results}/run.log) ...")

    cmd = build_snakemake_cmd(
        snakefile=WORKFLOW, workdir=workdir, configfile=configfile,
        results_dir_rel=rel_results, cores=args.cores, extra=extra,
        config_overrides=overrides,
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
        _print_footer(rel_results, elapsed, run_dir)
    return rc


if __name__ == "__main__":
    sys.exit(main())
