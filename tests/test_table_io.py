import sys
import warnings
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))

import pytest
from table_io import (
    decode_table_bytes, read_csv_rows, load_assay_table, read_metadata,
)

# The four non-ASCII strings from a real Excel-for-Mac assay table.
REAL_TEXT = (
    "assay,fwd,rev,probe,notes\n"
    "A1,ACGT,TTGA,,14–15 cycles; Fernández-No; 60 °C; Caamaño-Antelo\n"
)


@pytest.mark.parametrize("encoding", ["mac_roman", "cp1252"])
def test_decode_excel_codepages(encoding):
    data = REAL_TEXT.encode(encoding)
    with pytest.warns(UserWarning, match=encoding):
        text, enc = decode_table_bytes(data)
    assert enc == encoding
    assert text == REAL_TEXT


@pytest.mark.parametrize("encoding", ["utf-8", "utf-8-sig", "utf-16"])
def test_decode_unicode_without_warning(encoding):
    data = REAL_TEXT.encode(encoding)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        text, _ = decode_table_bytes(data)
    assert text == REAL_TEXT  # BOM stripped


def test_decode_cp1252_word_heuristic_typical_windows_text():
    # typical Windows Excel content: accented names + smart quotes
    s = "assay,notes\nA1,“José Müller” – résumé\n"
    with pytest.warns(UserWarning):
        text, enc = decode_table_bytes(s.encode("cp1252"))
    assert enc == "cp1252" and text == s


@pytest.mark.parametrize("delim", [",", ";", "\t"])
def test_read_csv_rows_delimiters(tmp_path, delim):
    p = tmp_path / "t.csv"
    p.write_text(delim.join(["assay", "fwd", "rev", "probe"]) + "\n"
                 + delim.join(["A1", "ACGT", "TTGA", ""]) + "\n")
    rows = read_csv_rows(p)
    assert rows == [{"assay": "A1", "fwd": "ACGT", "rev": "TTGA", "probe": ""}]


def test_read_csv_rows_strips_whitespace_and_nbsp_and_crlf(tmp_path):
    p = tmp_path / "t.csv"
    p.write_bytes("assay , fwd,rev,probe\r\nA1 , ACGT ,TTGA,\r\n".encode("utf-8"))
    assert read_csv_rows(p) == [{"assay": "A1", "fwd": "ACGT", "rev": "TTGA", "probe": ""}]


def test_read_csv_rows_quoted_commas(tmp_path):
    p = tmp_path / "t.csv"
    p.write_text('assay,fwd,rev,probe,notes\nA1,ACGT,TTGA,,"a, b, c"\n')
    assert read_csv_rows(p)[0]["notes"] == "a, b, c"


def _table(tmp_path, body, header="assay,fwd,rev,probe"):
    p = tmp_path / "assay_table.csv"
    p.write_text(header + "\n" + body)
    return p


def test_load_assay_table_valid_with_mods_and_degenerate(tmp_path):
    p = _table(tmp_path, "A1,ACGTRYN,ttga,/56-FAM/ACG/ZEN/TT[BHQ1]\nA2,ACGT,TTGA,\n")
    rows = load_assay_table(p)
    assert [r["assay"] for r in rows] == ["A1", "A2"]


def test_load_assay_table_mac_roman_file(tmp_path):
    p = tmp_path / "assay_table.csv"
    p.write_bytes(REAL_TEXT.encode("mac_roman"))
    with pytest.warns(UserWarning):
        rows = load_assay_table(p)
    assert rows[0]["notes"].startswith("14–15")


def test_load_assay_table_missing_column(tmp_path):
    p = _table(tmp_path, "A1,ACGT,TTGA\n", header="assay,fwd,rev")
    with pytest.raises(ValueError, match="missing required columns.*probe"):
        load_assay_table(p)


@pytest.mark.parametrize("body,match", [
    (",ACGT,TTGA,\n", "empty assay name"),
    ("A1,ACGT,TTGA,\nA1,ACGT,TTGA,\n", "duplicate assay name 'A1'"),
    ("Bac 16S,ACGT,TTGA,\n", "'Bac 16S'.*whitespace"),
    ("A1,,TTGA,\n", "'A1'.*fwd.*empty"),
    ("A1,ACGT,TT–GA,\n", "'A1'.*rev.*invalid"),
    ("A1,ACGT,TTGA,ACG T\n", "'A1'.*probe.*invalid"),
])
def test_load_assay_table_validation_errors(tmp_path, body, match):
    p = _table(tmp_path, body)
    with pytest.raises(ValueError, match=match):
        load_assay_table(p)


def test_read_metadata_mac_roman_semicolon(tmp_path):
    p = tmp_path / "metadata.csv"
    p.write_bytes("accession;organism_name\nGCF_1.1;Vibrio caamaño\n".encode("mac_roman"))
    with pytest.warns(UserWarning):
        df = read_metadata(p)
    assert list(df.columns) == ["accession", "organism_name"]
    assert df.loc[0, "organism_name"] == "Vibrio caamaño"
