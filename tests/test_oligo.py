import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))

import pytest
from oligo import parse_oligo, OligoError

TOMMY_BASES = "TGGAATCGTTTGACTGCATTT"


@pytest.mark.parametrize("raw", [
    "TG+GAATC+GTTT+GACTGCATTT",      # IDT '+' notation
    "TG[G]AATC[G]TTT[G]ACTGCATTT",   # bracket notation
    "TG[+G]AATC[+G]TTT[+G]ACTGCATTT",
    "tg+gaatc+gttt+gactgcattt",      # lowercase
])
def test_lna_notations_parse_identically(raw):
    o = parse_oligo(raw)
    assert o.bases == TOMMY_BASES
    assert o.lna == (2, 7, 11)


def test_plain_and_degenerate_sequence():
    o = parse_oligo("acgtRYSWKMBDHVN")
    assert o.bases == "ACGTRYSWKMBDHVN" and o.lna == () and o.notes == ()


def test_empty_is_empty():
    assert parse_oligo("").bases == ""


@pytest.mark.parametrize("raw,bases", [
    ("/56-FAM/ACG/ZEN/TT/3IABkFQ/", "ACGTT"),      # existing Vpop notation
    ("/5HEX/ACGT/3BHQ_1/", "ACGT"),
    ("ACGT[AmMC6]ACGT", "ACGTACGT"),                # multi-char bracket: non-base
    ("AC/iTAO/GT", "ACGT"),
    ("AC/iSpC3/GT", "ACGT"),
    ("ACG/iMe-dC/GT", "ACGCGT"),                    # base-bearing internal mods
    ("ACGT/ideoxyI/ACGT", "ACGTNACGT"),
    ("ACGTIACGT", "ACGTNACGT"),                     # bare inosine
    ("AC/ideoxyU/GT", "ACTGT"),
    ("ACGU", "ACGT"),
    ("AC/iBiodT/GT", "ACTGT"),
    ("AC/iFluorT/GT", "ACTGT"),
    ("A*C*G*T", "ACGT"),                            # phosphorothioate
])
def test_modifications(raw, bases):
    assert parse_oligo(raw).bases == bases


def test_notes_record_changes():
    o = parse_oligo("/56-FAM/A+CG/iMe-dC/T*/3IABkFQ/")
    text = " | ".join(o.notes)
    for frag in ("/56-FAM/", "/iMe-dC/", "/3IABkFQ/", "phosphorothioate"):
        assert frag in text


@pytest.mark.parametrize("raw,match", [
    ("AC/iMyDye/GT", r"unrecognized internal modification /iMyDye/"),
    ("AC/FOO/GT", r"unrecognized modification /FOO/"),
    ("AC/56-FAM/GT", r"5' modification /56-FAM/ must come first"),
    ("/3IABkFQ/ACGT", r"3' modification /3IABkFQ/ must come last"),
    ("ACG/ZEN", r"unclosed '/'"),
    ("ACG[BHQ1", r"unclosed '\['"),
    ("AC+NGT", r"'\+' must be followed by A, C, G or T"),
    ("AC+", r"'\+' must be followed by A, C, G or T"),
    ("AC[R]GT", r"\[R\].*LNA.*A, C, G or T"),
    ("AC[rA]GT", r"RNA.*not supported"),
    ("AC[mA]GT", r"2'-O-methyl.*not supported"),
    ("ACGrUACG", r"RNA.*not supported"),            # IDT rN in mixed case
    ("ACGmAACG", r"2'-O-methyl.*not supported"),
    ("ACG TT", r"invalid character"),
    ("AC–GT", r"invalid character"),
    ("#NAME?", r"Excel"),
    ("=GAATC", r"Excel"),
])
def test_errors(raw, match):
    with pytest.raises(OligoError, match=match):
        parse_oligo(raw)


def test_all_lowercase_iupac_r_m_are_degenerate_not_rna():
    # 'r'/'m' are IUPAC R/M; only flagged as RNA/2'-OMe in mixed-case sequences
    assert parse_oligo("acgrmacg").bases == "ACGRMACG"
    assert parse_oligo("ACGRMACG").bases == "ACGRMACG"


def test_oligo_error_is_value_error():
    assert issubclass(OligoError, ValueError)
