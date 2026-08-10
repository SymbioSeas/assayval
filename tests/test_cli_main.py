from pathlib import Path
from primeval import cli


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
        "assembly_dir: assemblies\nresults_dir: primeval_benchmarking_results\n"
    )
    captured = {}

    def fake_runner(cmd, logf):
        captured["cmd"] = cmd
        return 0

    rc = cli.main(["--run-name", "Vpop", "--directory", str(wd)], runner=fake_runner)
    assert rc == 0
    dirs = list((wd / "primeval_benchmarking_results").glob("Vpop_*"))
    assert len(dirs) == 1
    assert not (wd / "results").exists()
    assert any(a.startswith("results_dir=primeval_benchmarking_results/Vpop_")
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
