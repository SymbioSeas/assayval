"""LNA-aware scoring and per-assay thresholds in the detection engine."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))

import pandas as pd
import pytest
from detect import (
    lna_positions_match, filter_primer_hits, check_probe_in_amplicons, revcomp,
    run_detection_full,
)

PRIMER = "AGCCGAGCGTTACCAGC"            # 17 nt; index 5 is 'A'
PROBE = "ACGGGACAAAAAGGATGGCGAGTAC"     # 25 nt; index 4 is 'G'
REV = "CGAACGCAATGATTCTCTGAGC"


def _hit(sseq, sstart=100, send=116, qseqid="T_fwd", length=17):
    return dict(qseqid=qseqid, sseqid="c1", length=length, mismatch=0, gapopen=0,
                qstart=1, qend=length, sstart=sstart, send=send, sseq=sseq)


def _contig_for(hit, tail=30):
    sstart, send = int(hit['sstart']), int(hit['send'])
    core = hit['sseq'] if sstart <= send else revcomp(hit['sseq'])
    return "A" * (min(sstart, send) - 1) + core + "A" * tail


def _mutate(seq, i, base="T"):
    return seq[:i] + base + seq[i + 1:]


# --- lna_positions_match ------------------------------------------------------

def test_lna_positions_match():
    assert lna_positions_match(PRIMER, PRIMER, (5,))
    assert not lna_positions_match(PRIMER, _mutate(PRIMER, 5), (5,))
    assert lna_positions_match(PRIMER, _mutate(PRIMER, 6), (5,))   # mismatch elsewhere
    assert not lna_positions_match(PRIMER, _mutate(PRIMER, 5, "-"), (5,))  # contig edge
    assert lna_positions_match(PRIMER, PRIMER, ())


# --- primers ------------------------------------------------------------------

@pytest.mark.parametrize("strand", ["+", "-"])
def test_primer_mismatch_at_lna_rejected_when_exact(strand):
    sseq = _mutate(PRIMER, 5)                       # 1 mismatch, at the LNA
    hit = _hit(sseq) if strand == "+" else _hit(sseq, sstart=116, send=100)
    contigs = {"c1": _contig_for(hit)}
    kw = dict(primer_seq=PRIMER, max_mismatch=2, prime3_exact=1, contigs=contigs,
              lna_positions=(5,))
    assert len(filter_primer_hits(pd.DataFrame([hit]), **kw, lna_mismatch="exact")) == 0
    assert len(filter_primer_hits(pd.DataFrame([hit]), **kw, lna_mismatch="count")) == 1


def test_primer_mismatch_away_from_lna_still_tolerated():
    hit = _hit(_mutate(PRIMER, 8))
    out = filter_primer_hits(pd.DataFrame([hit]), PRIMER, 2, 1, {"c1": _contig_for(hit)},
                             lna_positions=(5,), lna_mismatch="exact")
    assert len(out) == 1


# --- probe --------------------------------------------------------------------

def _amp():
    return [{'contig_id': 'c1', 'amplicon_start': 100, 'amplicon_end': 300,
             'amplicon_size': 201, 'fwd_mismatches': 0, 'rev_mismatches': 0,
             'probe_found': False, 'probe_mismatches': None, 'probe_strand': None,
             'amplicon_sequence': '', 'strand': '+'}]


@pytest.mark.parametrize("policy,found", [("exact", False), ("count", True)])
def test_probe_mismatch_at_lna(policy, found):
    hit = _hit(_mutate(PROBE, 4), sstart=174, send=150, qseqid="T_probe", length=25)
    out = check_probe_in_amplicons(_amp(), pd.DataFrame([hit]), PROBE, 1,
                                   {"c1": _contig_for(hit)}, lna_positions=(4,),
                                   lna_mismatch=policy)
    assert out[0]['probe_found'] is found


# --- run_detection_full: notation, policy, per-assay thresholds --------------

def _setup(tmp_path, probe_cell, extra_cols="", extra_vals="", probe_site=PROBE):
    seq = ("T" * 99 + PRIMER + "T" * 33 + probe_site + "T" * 104
           + revcomp(REV) + "T" * 100)
    (tmp_path / "G.fna").write_text(f">c1\n{seq}\n")
    (tmp_path / "blast.tsv").write_text(
        "\t".join(["A_fwd", "c1", "17", "0", "1", "17", "100", "116"]) + "\n"
        + "\t".join(["A_rev", "c1", "22", "0", "1", "22", "300", "279"]) + "\n"
        + "\t".join(["A_probe", "c1", "25", "0", "1", "25", "150", "174"]) + "\n")
    (tmp_path / "assay_table.csv").write_text(
        f"assay,fwd,rev,probe{extra_cols}\nA,{PRIMER},{REV},{probe_cell}{extra_vals}\n")


def _call(tmp_path, **kw):
    args = dict(max_primer_mismatches=2, prime3_exact_nt=1, max_probe_mismatches=1,
                max_amplicon_size=500, store_amplicon_sequences=False)
    args.update(kw)
    det, _ = run_detection_full(str(tmp_path / "blast.tsv"),
                                str(tmp_path / "assay_table.csv"),
                                str(tmp_path / "G.fna"), **args)
    return det.iloc[0]['detection_call']


@pytest.mark.parametrize("probe_cell", ["ACGG+GACAAAAAGGATGGCGAGTAC",
                                        "ACGG[G]ACAAAAAGGATGGCGAGTAC"])
def test_lna_notation_end_to_end(tmp_path, probe_cell):
    # perfect probe site: both notations detect (bracket form used to lose a base)
    _setup(tmp_path, probe_cell)
    assert _call(tmp_path) == "Detected"
    # 1 mismatch AT the LNA: rejected under exact, tolerated under count
    _setup(tmp_path, probe_cell, probe_site=_mutate(PROBE, 4))
    assert _call(tmp_path) == "Primer Only"
    assert _call(tmp_path, lna_mismatch="count") == "Detected"


def test_per_assay_probe_threshold_overrides_config(tmp_path):
    # 1 probe mismatch (no LNA): config allows 1, the assay column demands 0
    _setup(tmp_path, PROBE, ",max_probe_mismatches", ",0",
           probe_site=_mutate(PROBE, 12))
    assert _call(tmp_path, max_probe_mismatches=1) == "Primer Only"
    _setup(tmp_path, PROBE, ",max_probe_mismatches", ",",   # blank -> config
           probe_site=_mutate(PROBE, 12))
    assert _call(tmp_path, max_probe_mismatches=1) == "Detected"


def test_per_assay_amplicon_size_overrides_config(tmp_path):
    _setup(tmp_path, PROBE, ",max_amplicon_size", ",150")    # amplicon is 201 bp
    assert _call(tmp_path, max_amplicon_size=500) == "Not Detected"
