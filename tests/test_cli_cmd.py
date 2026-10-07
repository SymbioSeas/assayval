from assayval.cli import build_parser, build_snakemake_cmd, filter_terminal_line


def test_parser_defaults():
    args = build_parser().parse_args([])
    assert args.run_name == "results"
    assert args.directory == "."
    assert args.cores == 8
    assert args.force is False
    assert args.configfile is None


def test_parser_passthrough_after_ddash():
    args = build_parser().parse_args(["--run-name", "Vpop", "--", "-n", "--rerun-incomplete"])
    assert args.run_name == "Vpop"
    # REMAINDER keeps the leading "--"; main() strips it
    assert args.snakemake_args[-2:] == ["-n", "--rerun-incomplete"]


def test_build_cmd_core_flags():
    cmd = build_snakemake_cmd(
        snakefile="/repo/workflow/Snakefile", workdir="/wd",
        configfile="/wd/config.yaml", results_dir_rel="results/Vpop_2026-07-07",
        cores=8,
    )
    assert cmd[0] == "snakemake"
    assert "--snakefile" in cmd and "/repo/workflow/Snakefile" in cmd
    assert "--directory" in cmd and "/wd" in cmd
    assert "--configfile" in cmd and "/wd/config.yaml" in cmd
    assert "--config" in cmd and "results_dir=results/Vpop_2026-07-07" in cmd
    assert cmd[cmd.index("--cores") + 1] == "8"
    assert cmd[cmd.index("--quiet") + 1] == "rules"


def test_build_cmd_extra_appended():
    cmd = build_snakemake_cmd(
        snakefile="s", workdir="w", configfile="c",
        results_dir_rel="results/x", cores=4, extra=["-n"],
    )
    assert cmd[-1] == "-n"


def test_filter_progress_milestones():
    state = {}
    assert filter_terminal_line("5 of 100 steps (5%) done", state)
    assert filter_terminal_line("8 of 100 steps (8%) done", state) == []
    assert filter_terminal_line("12 of 100 steps (12%) done", state)


def test_filter_alerts_pass_through():
    state = {}
    assert filter_terminal_line("Error in rule run_blast:", state)
    assert filter_terminal_line("WARNING: something", state)
    assert filter_terminal_line("Traceback (most recent call last):", state)


def test_filter_python_warnings_cleaned_and_shown_once():
    """Every per-assembly job re-reads the assay table, so a table warning
    arrives once per job; show it once, without the source-path prefix."""
    state = {}
    line = ("/x/workflow/scripts/table_io.py:102: UserWarning: assay_table.csv "
            "is not UTF-8; read it as Mac Roman\n")
    assert filter_terminal_line(line, state) == [
        "  warning: assay_table.csv is not UTF-8; read it as Mac Roman"]
    assert filter_terminal_line(line, state) == []
    # errors are never de-duplicated
    assert filter_terminal_line("Error in rule run_detection:", state)
    assert filter_terminal_line("Error in rule run_detection:", state)


def test_filter_suppresses_noise():
    state = {}
    assert filter_terminal_line("Building DAG of jobs...", state) == []
    assert filter_terminal_line("Select jobs to execute...", state) == []
    assert filter_terminal_line("Finished job 42.", state) == []


# --- re-score mode / config overrides ------------------------------------------

from pathlib import Path

import pytest

from assayval.cli import parse_overrides, resolve_rescore_blast_dir


def test_parser_rescore_defaults():
    args = build_parser().parse_args([])
    assert args.rescore_from is None
    assert args.overrides == []


def test_parser_collects_repeated_set_flags():
    args = build_parser().parse_args(
        ["--set", "max_primer_mismatches=1", "--set", "prime3_exact_nt=2"]
    )
    assert args.overrides == ["max_primer_mismatches=1", "prime3_exact_nt=2"]


def test_parse_overrides_passes_valid_tokens_through():
    tokens = ["max_primer_mismatches=1", "keep_blast=true"]
    assert parse_overrides(tokens) == tokens


@pytest.mark.parametrize("bad", ["max_primer_mismatches", "=2", ""])
def test_parse_overrides_rejects_malformed(bad):
    with pytest.raises(ValueError):
        parse_overrides([bad])


def test_build_cmd_overrides_share_the_single_config_flag():
    """A second --config would replace the first (nargs='+'), silently dropping
    results_dir, so overrides must be appended to the same flag."""
    cmd = build_snakemake_cmd(
        snakefile="s", workdir="w", configfile="c",
        results_dir_rel="results/x", cores=4,
        config_overrides=["blast_dir=results/prev/blast", "max_primer_mismatches=1"],
    )
    assert cmd.count("--config") == 1
    i = cmd.index("--config")
    assert cmd[i + 1:i + 4] == [
        "results_dir=results/x",
        "blast_dir=results/prev/blast",
        "max_primer_mismatches=1",
    ]
    # the override block must end before the next flag
    assert cmd[i + 4] == "--cores"


def _cache(tmp_path, name="prev_2026-08-19", n=2):
    blast = tmp_path / "results" / name / "blast"
    blast.mkdir(parents=True)
    for i in range(n):
        (blast / f"GCF_{i}.tsv.gz").write_bytes(b"")
    return blast


def test_resolve_rescore_accepts_run_directory(tmp_path):
    _cache(tmp_path)
    rel, n = resolve_rescore_blast_dir(tmp_path, "results/prev_2026-08-19")
    assert rel == "results/prev_2026-08-19/blast"
    assert n == 2


def test_resolve_rescore_accepts_blast_directory_directly(tmp_path):
    _cache(tmp_path)
    rel, n = resolve_rescore_blast_dir(tmp_path, "results/prev_2026-08-19/blast")
    assert rel == "results/prev_2026-08-19/blast"
    assert n == 2


def test_resolve_rescore_cache_outside_workdir_stays_absolute(tmp_path):
    other = tmp_path / "elsewhere"
    blast = other / "run" / "blast"
    blast.mkdir(parents=True)
    (blast / "GCF_0.tsv.gz").write_bytes(b"")
    wd = tmp_path / "analysis"
    wd.mkdir()
    rel, n = resolve_rescore_blast_dir(wd, str(blast))
    assert Path(rel).is_absolute()
    assert n == 1


def test_resolve_rescore_rejects_missing_path(tmp_path):
    with pytest.raises(ValueError, match="does not exist"):
        resolve_rescore_blast_dir(tmp_path, "results/nope")


def test_resolve_rescore_rejects_run_without_blast_dir(tmp_path):
    (tmp_path / "results" / "prev").mkdir(parents=True)
    with pytest.raises(ValueError, match="keep_blast: true"):
        resolve_rescore_blast_dir(tmp_path, "results/prev")


def test_resolve_rescore_rejects_empty_cache(tmp_path):
    (tmp_path / "results" / "prev" / "blast").mkdir(parents=True)
    with pytest.raises(ValueError, match="no cached BLAST output"):
        resolve_rescore_blast_dir(tmp_path, "results/prev")


def test_version_string_reports_package_and_commit():
    """Run provenance depends on this: a release number alone does not identify
    a working tree, so a git checkout should also surface its commit."""
    from assayval.cli import version_string
    v = version_string()
    assert isinstance(v, str) and v
    assert not v.startswith("(")


def test_parser_exposes_version_flag(capsys):
    import pytest as _pytest
    with _pytest.raises(SystemExit) as e:
        build_parser().parse_args(["--version"])
    assert e.value.code == 0
    assert "assay-val" in capsys.readouterr().out


# --- --set overrides reflected in the run header --------------------------------

def test_override_value_last_flag_wins():
    """Snakemake resolves repeated --config keys with the last value, so the CLI
    must report the same one or the header contradicts the run."""
    from assayval.cli import override_value
    ov = ["assembly_dir=first", "max_primer_mismatches=1", "assembly_dir=second"]
    assert override_value(ov, "assembly_dir") == "second"
    assert override_value(ov, "max_primer_mismatches") == "1"
    assert override_value(ov, "absent") is None


def test_count_assemblies_honours_assembly_dir_override(tmp_path):
    from assayval.cli import _count_assemblies
    (tmp_path / "assemblies").mkdir()
    for i in range(7):
        (tmp_path / "assemblies" / f"a{i}.fna").write_text(">x\nACGT\n")
    (tmp_path / "subset").mkdir()
    for i in range(2):
        (tmp_path / "subset" / f"a{i}.fna").write_text(">x\nACGT\n")
    cfg = tmp_path / "config.yaml"
    cfg.write_text("assembly_dir: assemblies\n")
    assert _count_assemblies(cfg, tmp_path) == 7
    assert _count_assemblies(cfg, tmp_path, "subset") == 2
