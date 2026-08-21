from pathlib import Path
from assayval import cli


def _make_workdir(tmp_path):
    (tmp_path / "assemblies").mkdir()
    (tmp_path / "config.yaml").write_text("assembly_dir: assemblies\nresults_dir: results\n")
    return tmp_path


def test_main_creates_run_dir_and_invokes_runner(tmp_path, capsys):
    wd = _make_workdir(tmp_path)
    captured = {}

    def fake_runner(cmd, logf):
        captured["cmd"] = cmd
        captured["logf"] = Path(logf)
        return 0

    rc = cli.main(["--run-name", "Vpop", "--directory", str(wd)], runner=fake_runner)
    assert rc == 0
    dirs = list((wd / "results").glob("Vpop_*"))
    assert len(dirs) == 1
    assert any(a.startswith("results_dir=results/Vpop_") for a in captured["cmd"])
    assert captured["logf"].name == "run.log"
    out = capsys.readouterr().out
    assert "run:" in out and "results/Vpop_" in out


def test_main_missing_config_errors(tmp_path, capsys):
    rc = cli.main(["--directory", str(tmp_path)], runner=lambda c, l: 0)
    assert rc == 2
    assert "config file not found" in capsys.readouterr().err


def test_main_reports_failure(tmp_path, capsys):
    wd = _make_workdir(tmp_path)
    rc = cli.main(["--directory", str(wd)], runner=lambda c, l: 1)
    assert rc == 1
    assert "run failed" in capsys.readouterr().err


def test_main_dry_run_does_not_create_dir(tmp_path):
    wd = _make_workdir(tmp_path)
    captured = {}

    def r(cmd, logf):
        captured["logf"] = logf
        return 0

    rc = cli.main(["--run-name", "Vpop", "--directory", str(wd), "--", "-n"], runner=r)
    assert rc == 0
    # dry run must not create a dated run directory or a run.log
    assert not list((wd / "results").glob("Vpop_*")) if (wd / "results").exists() else True
    assert captured["logf"] is None


def test_main_respects_results_dir_from_config(tmp_path):
    """config.yaml's results_dir must be used as the results root (it was
    previously ignored: the CLI hardcoded <workdir>/results)."""
    wd = tmp_path
    (wd / "assemblies").mkdir()
    (wd / "config.yaml").write_text(
        "assembly_dir: assemblies\nresults_dir: assayval_benchmarking_results\n"
    )
    captured = {}

    def fake_runner(cmd, logf):
        captured["cmd"] = cmd
        return 0

    rc = cli.main(["--run-name", "Vpop", "--directory", str(wd)], runner=fake_runner)
    assert rc == 0
    dirs = list((wd / "assayval_benchmarking_results").glob("Vpop_*"))
    assert len(dirs) == 1
    assert not (wd / "results").exists()
    assert any(a.startswith("results_dir=assayval_benchmarking_results/Vpop_")
               for a in captured["cmd"])


def test_main_results_dir_defaults_to_results(tmp_path):
    """A config without results_dir falls back to <workdir>/results."""
    wd = tmp_path
    (wd / "assemblies").mkdir()
    (wd / "config.yaml").write_text("assembly_dir: assemblies\n")

    rc = cli.main(["--run-name", "Vpop", "--directory", str(wd)],
                  runner=lambda c, l: 0)
    assert rc == 0
    assert len(list((wd / "results").glob("Vpop_*"))) == 1


# --- re-score mode -------------------------------------------------------------

def _make_cache(wd, name="prev_2026-08-19", n=2):
    blast = wd / "results" / name / "blast"
    blast.mkdir(parents=True)
    for i in range(n):
        (blast / f"GCF_{i}.tsv.gz").write_bytes(b"")
    return blast


def test_main_rescore_sets_blast_dir_and_thresholds(tmp_path, capsys):
    wd = _make_workdir(tmp_path)
    _make_cache(wd)
    captured = {}

    def fake_runner(cmd, logf):
        captured["cmd"] = cmd
        return 0

    rc = cli.main([
        "--run-name", "sweep_mm1", "--directory", str(wd),
        "--rescore-from", "results/prev_2026-08-19",
        "--set", "max_primer_mismatches=1",
    ], runner=fake_runner)
    assert rc == 0
    cmd = captured["cmd"]
    assert "blast_dir=results/prev_2026-08-19/blast" in cmd
    assert "max_primer_mismatches=1" in cmd
    assert cmd.count("--config") == 1
    out = capsys.readouterr().out
    assert "re-score" in out and "2 cached BLAST files" in out


def test_main_rescore_missing_cache_errors_before_creating_run_dir(tmp_path, capsys):
    wd = _make_workdir(tmp_path)
    (wd / "results" / "prev_2026-08-19").mkdir(parents=True)
    called = []

    rc = cli.main([
        "--run-name", "sweep", "--directory", str(wd),
        "--rescore-from", "results/prev_2026-08-19",
    ], runner=lambda c, l: called.append(c) or 0)
    assert rc == 2
    assert called == []
    assert "keep_blast: true" in capsys.readouterr().err
    # the failed invocation must not burn a dated run directory
    assert not list((wd / "results").glob("sweep_*"))


def test_main_rejects_malformed_set(tmp_path, capsys):
    wd = _make_workdir(tmp_path)
    rc = cli.main(["--directory", str(wd), "--set", "nonsense"], runner=lambda c, l: 0)
    assert rc == 2
    assert "KEY=VALUE" in capsys.readouterr().err


def test_main_header_reports_the_overridden_assembly_dir(tmp_path, capsys):
    """The calibration pilot runs with --set assembly_dir=<subset>; the header
    must count that directory, not the one named in the config file."""
    wd = _make_workdir(tmp_path)
    (wd / "assemblies" / "a.fna").write_text(">x\nACGT\n")
    (wd / "assemblies" / "b.fna").write_text(">x\nACGT\n")
    (wd / "subset").mkdir()
    (wd / "subset" / "a.fna").write_text(">x\nACGT\n")

    rc = cli.main(["--run-name", "pilot", "--directory", str(wd),
                   "--set", "assembly_dir=subset"], runner=lambda c, l: 0)
    assert rc == 0
    out = capsys.readouterr().out
    assert "subset/  (1 assemblies)" in out
