import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))

import pandas as pd
import pytest
from detect import (
    revcomp, find_valid_amplicons, extract_flanked_amplicon,
    run_detection_full, AMPLICON_RECORD_COLS,
)
from export_amplicons import (
    sanitize, organism_lookup, format_header, export_amplicon_fasta,
)

FWD = "AGCCGAGCGTTACCAGC"            # 17 nt
REV = "CGAACGCAATGATTCTCTGAGC"       # 22 nt
PROBE = "ACGGGACAAAAAGGATGGCGAGTAC"  # 25 nt


def _hit(qseqid, sstart, send, sseqid="c1", mismatch=0):
    return dict(qseqid=qseqid, sseqid=sseqid, length=17, gapopen=0, qstart=1,
                qend=17, sstart=sstart, send=send, mismatch=mismatch)


# --- strand + flanked extraction ---------------------------------------------

def test_find_valid_amplicons_records_strand():
    plus = find_valid_amplicons(pd.DataFrame([_hit("a_fwd", 100, 116)]),
                                pd.DataFrame([_hit("a_rev", 300, 279)]), 500)
    minus = find_valid_amplicons(pd.DataFrame([_hit("a_fwd", 2000, 1981)]),
                                 pd.DataFrame([_hit("a_rev", 1720, 1741)]), 500)
    assert plus[0]['strand'] == '+' and minus[0]['strand'] == '-'


def test_extract_flanked_plus_strand():
    contig = "A" * 10 + "CCGGTT" + "G" * 10          # amplicon at 11..16
    seq, up, down = extract_flanked_amplicon({"c": contig}, "c", 11, 16, "+", 3)
    assert seq == "aaa" + "CCGGTT" + "ggg"
    assert (up, down) == (3, 3)


def test_extract_flanked_minus_strand_is_oriented_fwd_to_rev():
    # genomic: left flank 'aaa', amplicon, right flank 'cct'; on the minus strand
    # the fwd primer sits at the HIGH coordinate, so the record must be the
    # reverse complement with the right-hand genomic flank upstream.
    contig = "TTTTAAA" + "CCGGTA" + "CCTTTTT"         # amplicon at 8..13
    seq, up, down = extract_flanked_amplicon({"c": contig}, "c", 8, 13, "-", 3)
    assert seq == revcomp("aaa" + "CCGGTA" + "cct")
    assert seq == "agg" + "TACCGG" + "ttt"
    assert (up, down) == (3, 3)


def test_extract_flanked_clips_at_contig_edges():
    contig = "AC" + "GGGG" + "T"                      # amplicon at 3..6
    seq, up, down = extract_flanked_amplicon({"c": contig}, "c", 3, 6, "+", 50)
    assert seq == "ac" + "GGGG" + "t"
    assert (up, down) == (2, 1)
    seq, up, down = extract_flanked_amplicon({"c": contig}, "c", 3, 6, "-", 50)
    assert (up, down) == (1, 2)


def test_extract_flanked_zero_flank():
    seq, up, down = extract_flanked_amplicon({"c": "AACCGGTT"}, "c", 3, 6, "+", 0)
    assert seq == "CCGG" and (up, down) == (0, 0)


# --- per-amplicon records from run_detection ---------------------------------

def _two_copy_genome(tmp_path, probe_in_first_only=True):
    """Two amplicon copies on one contig; copy 1 carries the probe site."""
    filler = "T" * 60
    copy1 = FWD + "T" * 20 + PROBE + "T" * 20 + revcomp(REV)
    copy2 = FWD + "T" * 65 + revcomp(REV)
    seq = filler + copy1 + filler + copy2 + filler
    (tmp_path / "GCF_9.1.fna").write_text(f">c1\n{seq}\n")
    s1 = 61
    e1 = s1 + len(copy1) - 1
    s2 = e1 + 61
    e2 = s2 + len(copy2) - 1
    p1 = s1 + len(FWD) + 20
    rows = [
        ["A_fwd", "c1", 17, 0, 1, 17, s1, s1 + 16],
        ["A_fwd", "c1", 17, 0, 1, 17, s2, s2 + 16],
        ["A_rev", "c1", 22, 0, 1, 22, e1, e1 - 21],
        ["A_rev", "c1", 22, 0, 1, 22, e2, e2 - 21],
        ["A_probe", "c1", 25, 0, 1, 25, p1, p1 + 24],
    ]
    (tmp_path / "blast.tsv").write_text(
        "".join("\t".join(map(str, r)) + "\n" for r in rows))
    (tmp_path / "assay_table.csv").write_text(
        f"assay,fwd,rev,probe\nA,{FWD},{REV},{PROBE}\n")
    return (s1, e1), (s2, e2)


def _run(tmp_path, flank=50, size=500):
    return run_detection_full(
        blast_tsv=str(tmp_path / "blast.tsv"),
        assay_table=str(tmp_path / "assay_table.csv"),
        fna_path=str(tmp_path / "GCF_9.1.fna"),
        max_primer_mismatches=0, prime3_exact_nt=1, max_probe_mismatches=0,
        max_amplicon_size=size, store_amplicon_sequences=False, flank_bp=flank)


def test_records_one_per_amplicon_in_detection_order(tmp_path):
    (s1, e1), (s2, e2) = _two_copy_genome(tmp_path)
    det, rec = _run(tmp_path, size=150)
    assert list(rec.columns) == AMPLICON_RECORD_COLS
    row = det.iloc[0]
    assert row['n_amplicons'] == 2
    assert list(rec['amplicon_index']) == [1, 2]
    assert list(rec['n_amplicons']) == [2, 2]
    # k-th record matches the k-th entry of the detection CSV's amplicon_starts
    assert [str(s) for s in rec['amplicon_start']] == row['amplicon_starts'].split(';')
    assert list(rec['probe_status']) == ['yes', 'no']
    assert set(rec['assembly_call']) == {'Detected'}
    r1 = rec.iloc[0]
    assert r1['sequence'].upper().startswith(("T" * 50) + FWD)
    assert r1['sequence'][:50].islower() and r1['sequence'][50:50 + 17] == FWD
    assert (r1['flank_up_bp'], r1['flank_down_bp']) == (50, 50)


def test_records_include_primer_only(tmp_path):
    _two_copy_genome(tmp_path)
    # drop the probe hit: probe-based assay, valid amplicons, no probe → Primer Only
    lines = (tmp_path / "blast.tsv").read_text().splitlines()
    (tmp_path / "blast.tsv").write_text("\n".join(l for l in lines if "probe" not in l) + "\n")
    det, rec = _run(tmp_path, size=150)
    assert det.iloc[0]['detection_call'] == 'Primer Only'
    assert len(rec) == 2 and set(rec['probe_status']) == {'no'}


def test_records_probe_free_assay_status_none(tmp_path):
    _two_copy_genome(tmp_path)
    (tmp_path / "assay_table.csv").write_text(f"assay,fwd,rev,probe\nA,{FWD},{REV},\n")
    _, rec = _run(tmp_path, size=150)
    assert set(rec['probe_status']) == {'none'}


def test_records_empty_when_no_amplicons(tmp_path):
    _two_copy_genome(tmp_path)
    det, rec = _run(tmp_path, size=50)   # both copies exceed max size
    assert det.iloc[0]['detection_call'] == 'Not Detected'
    assert rec.empty and list(rec.columns) == AMPLICON_RECORD_COLS


# --- export ------------------------------------------------------------------

def test_sanitize():
    assert sanitize("Bacillus cereus ATCC 14579") == "Bacillus_cereus_ATCC_14579"
    assert sanitize("E. coli O157:H7 (str. x/y)") == "E._coli_O157_H7_str._x_y"


def test_organism_lookup_fallbacks():
    meta = pd.DataFrame({"accession": ["G1", "G2"], "organism_name": ["Bacillus a", None],
                         "species": ["sp1", "sp2"]})
    assert organism_lookup(meta, ["species"]) == {"G1": "Bacillus a", "G2": "unknown"}
    no_org = meta.drop(columns="organism_name")
    assert organism_lookup(no_org, ["species"]) == {"G1": "sp1", "G2": "sp2"}
    assert organism_lookup(no_org.drop(columns="species"), []) == {"G1": "unknown", "G2": "unknown"}


def _rec(**kw):
    base = dict(accession="GCF_000007825.1", assay="Bac_16S", amplicon_index=1,
                n_amplicons=2, assembly_call="Detected", probe_status="yes",
                contig_id="NC_004722.1", amplicon_start=10234, amplicon_end=10311,
                strand="+", amplicon_size=78, flank_up_bp=50, flank_down_bp=50,
                fwd_mismatches=0, rev_mismatches=1, probe_mismatches=0,
                sequence="acgtACGTacgt")
    base.update(kw)
    return base


def test_format_header():
    fid, header = format_header(_rec(), "Bacillus cereus ATCC 14579")
    assert fid == "GCF_000007825.1_amp1of2_Bacillus_cereus_ATCC_14579"
    assert header == (
        ">GCF_000007825.1_amp1of2_Bacillus_cereus_ATCC_14579 "
        'organism="Bacillus cereus ATCC 14579" assay=Bac_16S assembly_call=Detected '
        "probe=yes contig=NC_004722.1:10234-10311(+) amplicon_bp=78 flank_bp=50,50 "
        "fwd_mm=0 rev_mm=1 probe_mm=0")


def test_format_header_missing_probe_mm():
    _, header = format_header(_rec(probe_status="none", probe_mismatches=None), "x")
    assert header.endswith("probe_mm=NA")


def test_export_writes_every_assay_including_empty(tmp_path):
    recs = pd.DataFrame([_rec(), _rec(amplicon_index=2, amplicon_start=20000,
                                      amplicon_end=20077, probe_status="no",
                                      probe_mismatches=None)],
                        columns=AMPLICON_RECORD_COLS)
    out = tmp_path / "amplicon_fasta"
    counts = export_amplicon_fasta(recs, ["Bac_16S", "Unused/Assay"],
                                   {"GCF_000007825.1": "Bacillus cereus"}, out)
    assert counts == {"Bac_16S": 2, "Unused/Assay": 0}
    fasta = (out / "Bac_16S_amplicons.fasta").read_text()
    assert fasta.count(">") == 2
    assert ">GCF_000007825.1_amp2of2_Bacillus_cereus " in fasta
    assert (out / "Unused_Assay_amplicons.fasta").read_text() == ""
    index = pd.read_csv(out / "amplicon_index.csv")
    assert list(index["fasta_id"]) == ["GCF_000007825.1_amp1of2_Bacillus_cereus",
                                       "GCF_000007825.1_amp2of2_Bacillus_cereus"]
    assert "sequence" not in index.columns
    assert set(index["fasta_file"]) == {"Bac_16S_amplicons.fasta"}


def test_export_rejects_filename_collision(tmp_path):
    empty = pd.DataFrame(columns=AMPLICON_RECORD_COLS)
    with pytest.raises(ValueError, match="same file name"):
        export_amplicon_fasta(empty, ["A/B", "A_B"], {}, tmp_path / "o")


def test_export_wraps_sequence_lines(tmp_path):
    recs = pd.DataFrame([_rec(sequence="A" * 150)], columns=AMPLICON_RECORD_COLS)
    export_amplicon_fasta(recs, ["Bac_16S"], {}, tmp_path / "o")
    lines = (tmp_path / "o" / "Bac_16S_amplicons.fasta").read_text().splitlines()
    assert [len(l) for l in lines[1:]] == [70, 70, 10]
