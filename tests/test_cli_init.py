import csv
from pathlib import Path

from assayval import cli

REPO = Path(__file__).resolve().parent.parent


def test_init_writes_templates_into_directory(tmp_path, capsys):
    rc = cli.main(["init", "--directory", str(tmp_path)])
    assert rc == 0
    cfg = tmp_path / "config.yaml"
    table = tmp_path / "assay_table.csv"
    assert cfg.exists() and table.exists()
    # config.yaml is the repo default, byte for byte (single source of truth)
    assert cfg.read_text() == (REPO / "config" / "config.yaml").read_text()
    out = capsys.readouterr().out
    assert "config.yaml" in out and "assay_table.csv" in out


def test_init_defaults_to_current_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    assert cli.main(["init"]) == 0
    assert (tmp_path / "config.yaml").exists()
    assert (tmp_path / "assay_table.csv").exists()


def test_init_assay_table_has_required_columns(tmp_path):
    cli.main(["init", "--directory", str(tmp_path)])
    with open(tmp_path / "assay_table.csv", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    assert {"assay", "fwd", "rev", "probe"} <= set(reader.fieldnames)
    assert rows, "template should include example rows showing the format"


def test_init_config_points_at_dropped_assay_table(tmp_path):
    import yaml
    cli.main(["init", "--directory", str(tmp_path)])
    cfg = yaml.safe_load((tmp_path / "config.yaml").read_text())
    assert cfg["assay_table"] == "assay_table.csv"


def test_init_does_not_overwrite_existing_files(tmp_path, capsys):
    (tmp_path / "config.yaml").write_text("mine: true\n")
    rc = cli.main(["init", "--directory", str(tmp_path)])
    assert rc == 1
    assert (tmp_path / "config.yaml").read_text() == "mine: true\n"
    # the missing file is still written
    assert (tmp_path / "assay_table.csv").exists()
    assert "--force" in capsys.readouterr().err


def test_init_force_overwrites(tmp_path):
    (tmp_path / "config.yaml").write_text("mine: true\n")
    rc = cli.main(["init", "--directory", str(tmp_path), "--force"])
    assert rc == 0
    assert "mine: true" not in (tmp_path / "config.yaml").read_text()


def test_init_creates_missing_directory(tmp_path):
    target = tmp_path / "new_project"
    assert cli.main(["init", "--directory", str(target)]) == 0
    assert (target / "config.yaml").exists()


def test_run_mode_still_parses_without_init(tmp_path, capsys):
    # the flat run interface is unchanged: no config -> the usual error,
    # which now points users at `assay-val init`
    rc = cli.main(["--directory", str(tmp_path)], runner=lambda c, l: 0)
    assert rc == 2
    assert "assay-val init" in capsys.readouterr().err


def test_init_assay_table_has_utf8_bom(tmp_path):
    """The BOM makes Excel open the file as "CSV UTF-8", so a plain Save keeps
    UTF-8 instead of falling back to Mac Roman / Windows-1252."""
    cli.main(["init", "--directory", str(tmp_path)])
    assert (tmp_path / "assay_table.csv").read_bytes().startswith(b"\xef\xbb\xbf")


def test_config_template_is_pure_ascii():
    """An ASCII config survives any editor's re-encoding unchanged."""
    data = (REPO / "config" / "config.yaml").read_bytes()
    assert all(b < 128 for b in data), "non-ASCII byte in config/config.yaml"


def test_template_lists_per_assay_threshold_columns(tmp_path):
    cli.main(["init", "--directory", str(tmp_path)])
    header = (tmp_path / "assay_table.csv").read_text(encoding="utf-8-sig").splitlines()[0]
    for col in ("max_primer_mismatches", "prime3_exact_nt",
                "max_probe_mismatches", "max_amplicon_size"):
        assert col in header.split(",")
