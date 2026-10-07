import sys
import pytest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))

from prepare_oligos import strip_idt_modifications, load_assay_table, write_oligo_fasta


def test_strip_zen():
    assert strip_idt_modifications("CGTATTACC/ZEN/GCGGCTGCTGGCAC") == "CGTATTACCGCGGCTGCTGGCAC"


def test_strip_no_modification():
    assert strip_idt_modifications("AGCCGAGCGTTACCAGC") == "AGCCGAGCGTTACCAGC"


def test_strip_multiple_modifications():
    assert strip_idt_modifications("/56-FAM/ACGT/3IABkFQ/") == "ACGT"


def test_strip_internal_zen():
    assert strip_idt_modifications("AGCGCACAT/ZEN/CAGAAGTCGGCCA") == "AGCGCACATCAGAAGTCGGCCA"


def test_strip_bracket_modification():
    assert strip_idt_modifications("ACGT[AmMC6]ACGT") == "ACGTACGT"


def test_load_assay_table(tmp_path):
    csv_path = tmp_path / "assays.csv"
    csv_path.write_text("assay,probe,fwd,rev\nMyAssay,ACGT,TTTT,GGGG\n")
    assays = load_assay_table(str(csv_path))
    assert len(assays) == 1
    assert assays[0]['assay'] == 'MyAssay'
    assert assays[0]['probe'] == 'ACGT'
    assert assays[0]['fwd'] == 'TTTT'


def test_load_assay_table_missing_column_raises(tmp_path):
    csv_path = tmp_path / "bad.csv"
    csv_path.write_text("assay,probe,fwd\nX,ACGT,TTTT\n")  # missing 'rev'
    with pytest.raises(ValueError, match="missing required columns"):
        load_assay_table(str(csv_path))


def test_load_assay_table_bom(tmp_path):
    csv_path = tmp_path / "bom.csv"
    # Write file with UTF-8 BOM (as Excel would produce)
    csv_path.write_bytes(b'\xef\xbb\xbfassay,probe,fwd,rev\nX,ACGT,TTTT,CCCC\n')
    assays = load_assay_table(str(csv_path))
    assert assays[0]['assay'] == 'X'  # No leading BOM character in key


def test_write_fasta_strips_modifications(tmp_path):
    assays = [{'assay': 'TestA', 'probe': 'ACGT/ZEN/ACGT', 'fwd': 'TTTT', 'rev': 'CCCC'}]
    fasta_out = tmp_path / "oligos.fasta"
    log_out = tmp_path / "oligos.log"
    write_oligo_fasta(assays, str(fasta_out), str(log_out))
    content = fasta_out.read_text()
    assert '>TestA_fwd\nTTTT\n' in content
    assert '>TestA_rev\nCCCC\n' in content
    assert '>TestA_probe\nACGTACGT\n' in content
    assert '/ZEN/' not in content


def test_write_fasta_log_records_stripped(tmp_path):
    assays = [{'assay': 'X', 'probe': 'A/ZEN/T', 'fwd': 'GG', 'rev': 'CC'}]
    fasta_out = tmp_path / "oligos.fasta"
    log_out = tmp_path / "oligos.log"
    write_oligo_fasta(assays, str(fasta_out), str(log_out))
    log_content = log_out.read_text()
    assert 'X_probe' in log_content
    assert '/ZEN/' in log_content


def test_write_fasta_correct_count(tmp_path):
    assays = [
        {'assay': 'A1', 'probe': 'ACGT', 'fwd': 'TTTT', 'rev': 'CCCC'},
        {'assay': 'A2', 'probe': 'GGGG', 'fwd': 'AAAA', 'rev': 'CCCC'},
    ]
    fasta_out = tmp_path / "oligos.fasta"
    log_out = tmp_path / "oligos.log"
    write_oligo_fasta(assays, str(fasta_out), str(log_out))
    headers = [l for l in fasta_out.read_text().splitlines() if l.startswith('>')]
    assert len(headers) == 6  # 2 assays × 3 oligos


def test_load_assay_table_excel_mac_roman(tmp_path):
    """Excel for Mac saves plain CSV as Mac Roman; non-ASCII free text in a
    notes column must not crash the run (it used to: UnicodeDecodeError)."""
    csv_path = tmp_path / "assay_table.csv"
    csv_path.write_bytes(
        "assay,fwd,rev,probe,notes\r\nX,TTTT,CCCC,,14–15 cycles at 60 °C\r\n"
        .encode("mac_roman"))
    with pytest.warns(UserWarning, match="Mac Roman"):
        assays = load_assay_table(str(csv_path))
    assert assays[0]['notes'] == "14–15 cycles at 60 °C"


def test_load_assay_table_trims_trailing_space_in_primer(tmp_path):
    csv_path = tmp_path / "assay_table.csv"
    csv_path.write_text("assay,fwd,rev,probe\nX,TTTT ,CCCC,\n")
    assert load_assay_table(str(csv_path))[0]['fwd'] == 'TTTT'
