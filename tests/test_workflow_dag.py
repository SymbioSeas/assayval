"""Snakemake dry-run checks of the workflow DAG (no BLAST is run)."""
import shutil
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SNAKEFILE = REPO / "workflow" / "Snakefile"

pytestmark = pytest.mark.skipif(shutil.which("snakemake") is None,
                                reason="snakemake not on PATH")


def _analysis_dir(tmp_path):
    (tmp_path / "assemblies").mkdir()
    for acc in ("GCF_1.1", "GCF_2.1"):
        (tmp_path / "assemblies" / f"{acc}.fna").write_text(">c1\nACGT\n")
    (tmp_path / "assemblies" / "metadata.csv").write_text(
        "accession,organism_name\nGCF_1.1,A a\nGCF_2.1,B b\n")
    (tmp_path / "assay_table.csv").write_text("assay,fwd,rev,probe\nA,ACGT,ACGT,\n")
    shutil.copy(REPO / "config" / "config.yaml", tmp_path / "config.yaml")
    return tmp_path


def _dry_run(wd, *config):
    proc = subprocess.run(
        ["snakemake", "-n", "--snakefile", str(SNAKEFILE), "--directory", str(wd),
         "--configfile", str(wd / "config.yaml"), "--cores", "1",
         "--config", "results_dir=results/t", "group_by=organism_name", *config],
        capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return proc.stdout + proc.stderr


def test_dag_includes_amplicon_export_by_default(tmp_path):
    out = _dry_run(_analysis_dir(tmp_path))
    assert "export_amplicon_fasta" in out


def test_dag_omits_amplicon_export_when_disabled(tmp_path):
    out = _dry_run(_analysis_dir(tmp_path), "amplicon_fasta=false")
    assert "export_amplicon_fasta" not in out
    assert "run_detection" in out
