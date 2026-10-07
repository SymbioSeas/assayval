import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "workflow" / "scripts"))

import pytest
import pandas as pd
from detect import (
    check_3prime_exact,
    filter_primer_hits,
    find_valid_amplicons,
    check_probe_in_amplicons,
    call_detection,
    load_blast_results,
    iupac_match,
    count_iupac_mismatches,
    count_degenerate,
    revcomp,
    load_contigs,
    reconstruct_full_hits,
    BLAST_COLS,
)


def make_hit(**kwargs) -> dict:
    """BLAST hit dict with sensible defaults (fwd primer, + strand, perfect match)."""
    defaults = {
        'qseqid': 'TestAssay_fwd', 'sseqid': 'contig1',
        'pident': 100.0, 'length': 17, 'mismatch': 0, 'gapopen': 0,
        'qstart': 1, 'qend': 17, 'sstart': 100, 'send': 116,
        'evalue': 0.001, 'bitscore': 32.0,
        'qseq': 'AGCCGAGCGTTACCAGC', 'sseq': 'AGCCGAGCGTTACCAGC',
    }
    defaults.update(kwargs)
    return defaults


def contig_embedding(hit, pad_char='A', tail=30):
    """Build a contig in which the hit's sseq occupies its subject coordinates,
    so reconstruction re-derives exactly that window."""
    sstart, send = int(hit['sstart']), int(hit['send'])
    lo = min(sstart, send)
    core = hit['sseq'] if sstart <= send else revcomp(hit['sseq'])
    return pad_char * (lo - 1) + core + pad_char * tail


# --- check_3prime_exact ---

def test_3prime_perfect_match():
    assert check_3prime_exact('AGCCGAGCGTTACCAGC', 'AGCCGAGCGTTACCAGC') is True


def test_3prime_last_nt_mismatch():
    assert check_3prime_exact('AGCCGAGCGTTACCAGC', 'AGCCGAGCGTTACCAGT') is False


def test_3prime_third_from_end_mismatch():
    assert check_3prime_exact('AGCCGAGCGTTACCAGC', 'AGCCGAGCGTTACCTGC') is False


def test_3prime_gap_fails():
    assert check_3prime_exact('AGCCGAGCGTTACCAGC', 'AGCCGAGCGTTACCA-C') is False


def test_3prime_mismatch_outside_window_passes():
    assert check_3prime_exact('AGCCGAGCGTTACCAGC', 'TGCCGAGCGTTACCAGC') is True


# --- filter_primer_hits ---

def test_filter_primer_perfect_forward():
    hit = make_hit(qend=17, sstart=100, send=116, mismatch=0)
    contigs = {'contig1': contig_embedding(hit)}
    hits = pd.DataFrame([hit])
    result = filter_primer_hits(hits, primer_seq='AGCCGAGCGTTACCAGC',
                                max_mismatch=2, prime3_exact=3, contigs=contigs)
    assert len(result) == 1


def test_filter_primer_too_many_mismatches():
    # sseq has 3 mismatches vs primer (positions 0, 5, 10 changed)
    hit = make_hit(qend=17, mismatch=3, sstart=100, send=116,
                   sseq='TGCCGTGCGTTGCCAGC')
    contigs = {'contig1': contig_embedding(hit)}
    result = filter_primer_hits(pd.DataFrame([hit]), primer_seq='AGCCGAGCGTTACCAGC',
                                max_mismatch=2, prime3_exact=3, contigs=contigs)
    assert len(result) == 0


def test_filter_recovers_5prime_trimmed_hit():
    """THE FIX: a hit BLAST trimmed at the 5' end (2 mismatches at primer
    positions 2-3, clean 3' end) is reconstructed and passes."""
    primer = 'GATTCGGTGAAGAAGAGATGATCTC'          # 25 nt (groEL_Valg fwd)
    window = 'GTATCGGTGAAGAAGAGATGATCTC'          # genome, mismatches at pos 2,3
    contigs = {'contig1': 'C' * 99 + window + 'C' * 50}
    hit = make_hit(qstart=4, qend=25, length=22, mismatch=0,
                   sstart=103, send=124, qseq=primer[3:], sseq=window[3:])
    result = filter_primer_hits(pd.DataFrame([hit]), primer_seq=primer,
                                max_mismatch=3, prime3_exact=2, contigs=contigs)
    assert len(result) == 1
    assert result.iloc[0]['mismatch'] == 2
    assert (result.iloc[0]['sstart'], result.iloc[0]['send']) == (100, 124)


def test_filter_rejects_reconstructed_3prime_mismatch():
    """Tall_dnaJ scenario: reconstruction recovers a 3'-trimmed hit, but the
    penultimate-base mismatch still fails the 3'-exact filter."""
    primer = 'AGCAGCTTATGACCAATACGCC'             # 22 nt
    window = 'AGCAGCTTATGACCAATACGGC'             # pos 21 mismatch
    contigs = {'contig1': 'A' * 99 + window + 'A' * 30}
    hit = make_hit(qstart=1, qend=18, length=18, mismatch=0,
                   sstart=100, send=117, qseq=primer[:18], sseq=window[:18])
    result = filter_primer_hits(pd.DataFrame([hit]), primer_seq=primer,
                                max_mismatch=3, prime3_exact=2, contigs=contigs)
    assert len(result) == 0


def test_filter_rejects_trimmed_hit_over_mismatch_budget():
    """A spurious core hit reconstructs to a window with many mismatches and
    is rejected by the mismatch budget (length 10 passes the pre-filter,
    min_len = 17 - 4*2 = 9, so the reconstruction path is exercised)."""
    primer = 'AGCCGAGCGTTACCAGC'                  # 17 nt
    window = 'TTTTTTTTTTTACCAGC'                  # first 9 primer bases mismatch
    contigs = {'contig1': 'G' * 99 + window + 'G' * 30}
    hit = make_hit(qstart=8, qend=17, length=10, mismatch=2,
                   sstart=107, send=116, qseq=primer[7:], sseq=window[7:])
    result = filter_primer_hits(pd.DataFrame([hit]), primer_seq=primer,
                                max_mismatch=2, prime3_exact=3, contigs=contigs)
    assert len(result) == 0


def test_filter_prefilter_drops_short_spurious_cores():
    """Hits shorter than primer_len - 4*(budget + degenerate) cannot pass the
    mismatch filter (blastn-short trimming implies >= k/4 mismatches per
    k-base trimmed terminus) and are dropped before reconstruction."""
    primer = 'AGCCGAGCGTTACCAGC'                  # 17 nt; min_len = 9
    hit = make_hit(qstart=11, qend=17, length=7, mismatch=0,
                   sstart=110, send=116, qseq=primer[10:], sseq=primer[10:])
    result = filter_primer_hits(pd.DataFrame([hit]), primer_seq=primer,
                                max_mismatch=2, prime3_exact=3,
                                contigs={'contig1': 'G' * 200})
    assert len(result) == 0


def test_count_degenerate():
    assert count_degenerate('ACGT') == 0
    assert count_degenerate('CAGGTTTGYTGCACGGCGAAGA') == 1
    assert count_degenerate('GATCGAAGTRCCRACACTMGGA') == 3


def test_filter_primer_both_strands_accepted():
    """filter_primer_hits accepts both + and - strand hits; strand pairing is
    done in find_valid_amplicons. The two sites live at different coordinates
    of one contig."""
    primer = 'AGCCGAGCGTTACCAGC'
    plus_hit = make_hit(qend=17, mismatch=0, sstart=100, send=116)
    minus_hit = make_hit(qend=17, mismatch=0, sstart=316, send=300)
    contig = 'A' * 99 + primer + 'A' * 183 + revcomp(primer) + 'A' * 30
    contigs = {'contig1': contig}
    result = filter_primer_hits(pd.DataFrame([plus_hit, minus_hit]), primer_seq=primer,
                                max_mismatch=2, prime3_exact=3, contigs=contigs)
    assert len(result) == 2


def test_filter_primer_minus_strand_accepted():
    hit = make_hit(qend=22, mismatch=0, sstart=300, send=279,
                   qseq='CGAACGCAATGATTCTCTGAGC',
                   sseq='CGAACGCAATGATTCTCTGAGC')
    contigs = {'contig1': contig_embedding(hit)}
    result = filter_primer_hits(pd.DataFrame([hit]), primer_seq='CGAACGCAATGATTCTCTGAGC',
                                max_mismatch=2, prime3_exact=3, contigs=contigs)
    assert len(result) == 1


def test_filter_primer_3prime_mismatch_rejected():
    hit = make_hit(qend=17, mismatch=1,
                   qseq='AGCCGAGCGTTACCAGT',
                   sseq='AGCCGAGCGTTACCAGC',
                   sstart=100, send=116)
    contigs = {'contig1': contig_embedding(hit)}
    result = filter_primer_hits(pd.DataFrame([hit]), primer_seq='AGCCGAGCGTTACCAGT',
                                max_mismatch=2, prime3_exact=3, contigs=contigs)
    assert len(result) == 0


# --- find_valid_amplicons ---

def test_find_amplicons_valid():
    fwd = pd.DataFrame([make_hit(sseqid='c1', sstart=100, send=116, mismatch=0)])
    rev = pd.DataFrame([make_hit(qseqid='TestAssay_rev', sseqid='c1',
                                  sstart=300, send=279, mismatch=0)])
    amps = find_valid_amplicons(fwd, rev, max_amplicon_size=500)
    assert len(amps) == 1
    assert amps[0]['amplicon_size'] == 201  # 300 - 100 + 1


def test_find_amplicons_exceeds_size():
    fwd = pd.DataFrame([make_hit(sseqid='c1', sstart=100, send=116, mismatch=0)])
    rev = pd.DataFrame([make_hit(qseqid='TestAssay_rev', sseqid='c1',
                                  sstart=700, send=679, mismatch=0)])
    amps = find_valid_amplicons(fwd, rev, max_amplicon_size=500)
    assert len(amps) == 0


def test_find_amplicons_different_contigs():
    fwd = pd.DataFrame([make_hit(sseqid='c1', sstart=100, send=116, mismatch=0)])
    rev = pd.DataFrame([make_hit(qseqid='TestAssay_rev', sseqid='c2',
                                  sstart=300, send=279, mismatch=0)])
    amps = find_valid_amplicons(fwd, rev, max_amplicon_size=500)
    assert len(amps) == 0


def test_find_amplicons_inverted_rejected():
    fwd = pd.DataFrame([make_hit(sseqid='c1', sstart=300, send=316, mismatch=0)])
    rev = pd.DataFrame([make_hit(qseqid='TestAssay_rev', sseqid='c1',
                                  sstart=150, send=129, mismatch=0)])
    amps = find_valid_amplicons(fwd, rev, max_amplicon_size=500)
    assert len(amps) == 0


def test_find_amplicons_multi():
    fwd = pd.DataFrame([make_hit(sseqid='c1', sstart=100, send=116, mismatch=0)])
    rev = pd.DataFrame([
        make_hit(qseqid='TestAssay_rev', sseqid='c1', sstart=300, send=279, mismatch=0),
        make_hit(qseqid='TestAssay_rev', sseqid='c1', sstart=450, send=429, mismatch=1,
                 qseq='CGAACGCAATGATTCTCTGAGC', sseq='CGAACGCAATGATTCTCTGAGC'),
    ])
    amps = find_valid_amplicons(fwd, rev, max_amplicon_size=500)
    assert len(amps) == 2


def test_find_amplicons_reverse_orientation():
    """Target gene on - strand: fwd hits - strand, rev hits + strand → valid amplicon."""
    # fwd on - strand: sstart=2000, send=1981 (sstart > send)
    fwd = pd.DataFrame([make_hit(sseqid='c1', sstart=2000, send=1981, mismatch=0,
                                  qseq='AGCCGAGCGTTACCAGC', sseq='AGCCGAGCGTTACCAGC')])
    # rev on + strand: sstart=1720, send=1741 (sstart < send)
    rev = pd.DataFrame([make_hit(qseqid='TestAssay_rev', sseqid='c1',
                                  sstart=1720, send=1741, mismatch=0, qend=22,
                                  qseq='CGAACGCAATGATTCTCTGAGC',
                                  sseq='CGAACGCAATGATTCTCTGAGC')])
    amps = find_valid_amplicons(fwd, rev, max_amplicon_size=500)
    assert len(amps) == 1
    assert amps[0]['amplicon_start'] == 1720  # min(2000, 1720)
    assert amps[0]['amplicon_end'] == 2000    # max(2000, 1720)
    assert amps[0]['amplicon_size'] == 281    # 2000 - 1720 + 1


def test_find_amplicons_same_strand_rejected():
    """Both primers on the same strand cannot form an amplicon."""
    fwd = pd.DataFrame([make_hit(sseqid='c1', sstart=100, send=116, mismatch=0)])
    rev = pd.DataFrame([make_hit(qseqid='TestAssay_rev', sseqid='c1',
                                  sstart=300, send=316, mismatch=0)])  # both + strand
    amps = find_valid_amplicons(fwd, rev, max_amplicon_size=500)
    assert len(amps) == 0


# --- check_probe_in_amplicons ---

def _make_amplicon(contig='c1', start=100, end=300, fwd_mm=0, rev_mm=0):
    return {
        'contig_id': contig, 'amplicon_start': start, 'amplicon_end': end,
        'amplicon_size': end - start + 1, 'fwd_mismatches': fwd_mm,
        'rev_mismatches': rev_mm, 'probe_found': False,
        'probe_mismatches': None, 'probe_strand': None, 'amplicon_sequence': '',
    }


_PROBE_SEQ = 'ACGGGACAAAAAGGATGGCGAGTAC'  # 25 nt


def test_probe_within_amplicon():
    amps = [_make_amplicon()]
    hit = make_hit(qseqid='TestAssay_probe', sseqid='c1',
                   qend=25, length=25, sstart=150, send=174, mismatch=0,
                   qseq=_PROBE_SEQ, sseq=_PROBE_SEQ)
    contigs = {'c1': contig_embedding(hit)}
    result = check_probe_in_amplicons(amps, pd.DataFrame([hit]), probe_seq=_PROBE_SEQ,
                                      max_probe_mismatches=1, contigs=contigs)
    assert result[0]['probe_found'] is True
    assert result[0]['probe_mismatches'] == 0
    assert result[0]['probe_strand'] == '+'


def test_probe_outside_amplicon():
    amps = [_make_amplicon()]
    hit = make_hit(qseqid='TestAssay_probe', sseqid='c1',
                   qend=25, length=25, sstart=500, send=524, mismatch=0,
                   qseq=_PROBE_SEQ, sseq=_PROBE_SEQ)
    contigs = {'c1': contig_embedding(hit)}
    result = check_probe_in_amplicons(amps, pd.DataFrame([hit]), probe_seq=_PROBE_SEQ,
                                      max_probe_mismatches=1, contigs=contigs)
    assert result[0]['probe_found'] is False


def test_probe_too_many_mismatches():
    amps = [_make_amplicon()]
    # 2 true mismatches in sseq (positions 23,24 changed) so IUPAC count == 2
    sseq_mm2 = _PROBE_SEQ[:-2] + 'TT'
    hit = make_hit(qseqid='TestAssay_probe', sseqid='c1',
                   qend=25, length=25, sstart=150, send=174, mismatch=2,
                   qseq=_PROBE_SEQ, sseq=sseq_mm2)
    contigs = {'c1': contig_embedding(hit)}
    result = check_probe_in_amplicons(amps, pd.DataFrame([hit]), probe_seq=_PROBE_SEQ,
                                      max_probe_mismatches=1, contigs=contigs)
    assert result[0]['probe_found'] is False


def test_probe_minus_strand_detected():
    amps = [_make_amplicon()]
    hit = make_hit(qseqid='TestAssay_probe', sseqid='c1',
                   qend=25, length=25, sstart=174, send=150, mismatch=0,
                   qseq=_PROBE_SEQ, sseq=_PROBE_SEQ)
    contigs = {'c1': contig_embedding(hit)}
    result = check_probe_in_amplicons(amps, pd.DataFrame([hit]), probe_seq=_PROBE_SEQ,
                                      max_probe_mismatches=1, contigs=contigs)
    assert result[0]['probe_found'] is True
    assert result[0]['probe_strand'] == '-'


def test_probe_5prime_trimmed_recovered():
    """Probe hits get the same end-trimming reconstruction as primers."""
    amps = [_make_amplicon()]
    window = 'T' + _PROBE_SEQ[1:]                 # mismatch at probe position 1
    hit = make_hit(qseqid='TestAssay_probe', sseqid='c1',
                   qstart=2, qend=25, length=24, mismatch=0,
                   sstart=151, send=174, qseq=_PROBE_SEQ[1:], sseq=window[1:])
    contigs = {'c1': 'G' * 149 + window + 'G' * 30}
    result = check_probe_in_amplicons(amps, pd.DataFrame([hit]), probe_seq=_PROBE_SEQ,
                                      max_probe_mismatches=1, contigs=contigs)
    assert result[0]['probe_found'] is True
    assert result[0]['probe_mismatches'] == 1


# --- call_detection ---

def test_call_detected():
    amps = [{'probe_found': True, 'amplicon_size': 142, 'contig_id': 'c1'}]
    call, n, multi, sizes, contigs = call_detection(amps)
    assert call == 'Detected'
    assert n == 1
    assert multi is False
    assert '142' in sizes


def test_call_primer_only():
    amps = [{'probe_found': False, 'amplicon_size': 142, 'contig_id': 'c1'}]
    call, n, multi, sizes, contigs = call_detection(amps)
    assert call == 'Primer Only'
    assert n == 1


def test_call_not_detected():
    call, n, multi, sizes, contigs = call_detection([])
    assert call == 'Not Detected'
    assert n == 0
    assert sizes == ''


def test_call_multi_amplicon():
    amps = [
        {'probe_found': True, 'amplicon_size': 142, 'contig_id': 'c1'},
        {'probe_found': True, 'amplicon_size': 198, 'contig_id': 'c1'},
    ]
    call, n, multi, sizes, contigs = call_detection(amps)
    assert call == 'Detected'
    assert n == 2
    assert multi is True
    assert '142' in sizes
    assert '198' in sizes


# --- load_blast_results ---

def test_load_blast_empty(tmp_path):
    f = tmp_path / "empty.tsv"
    f.write_text('')
    df = load_blast_results(str(f))
    assert df.empty
    assert list(df.columns) == BLAST_COLS


def test_load_blast_nonexistent_file(tmp_path):
    df = load_blast_results(str(tmp_path / "no_such_file.tsv"))
    assert df.empty
    assert list(df.columns) == BLAST_COLS


# --- run_detection integration ---

def test_run_detection_integration(tmp_path):
    """End-to-end test: one assay, one assembly, full detection."""
    from detect import run_detection

    fwd_seq = "AGCCGAGCGTTACCAGC"                  # 17 nt, at 100..116 (+)
    rev_seq = "CGAACGCAATGATTCTCTGAGC"             # 22 nt, revcomp at 279..300
    probe_seq = "ACGGGACAAAAAGGATGGCGAGTAC"        # 25 nt, at 150..174 (+)

    # Minimal assay table — VhPath only, no IDT modifications in these seqs
    assay_csv = tmp_path / "assay_table.csv"
    assay_csv.write_text(
        "assay,probe,fwd,rev\n"
        f"VhPath,{probe_seq},{fwd_seq},{rev_seq}\n"
    )

    # .fna genuinely containing the amplicon at positions 100-300:
    # 1..99 filler | fwd 100..116 | filler | probe 150..174 | filler |
    # revcomp(rev) 279..300 | filler to 400. Hit windows are re-extracted
    # from this sequence by reconstruct_full_hits, so it must be real.
    seq = ("T" * 99 + fwd_seq + "T" * 33 + probe_seq + "T" * 104
           + revcomp(rev_seq) + "T" * 100)
    assert len(seq) == 400
    fna_path = tmp_path / "GCF_000001.fna"
    fna_path.write_text(f">contig1\n{seq}\n")

    # BLAST TSV: fwd hit (+strand), rev hit (-strand), probe hit within amplicon
    # Columns are qseqid sseqid length gapopen qstart qend sstart send, matching
    # BLAST_COLS and the -outfmt in workflow/rules/blast.smk.
    blast_tsv = tmp_path / "blast.tsv"
    blast_tsv.write_text(
        "\t".join(["VhPath_fwd", "contig1", "17", "0", "1", "17", "100", "116"]) + "\n" +
        "\t".join(["VhPath_rev", "contig1", "22", "0", "1", "22", "300", "279"]) + "\n" +
        "\t".join(["VhPath_probe", "contig1", "25", "0", "1", "25", "150", "174"]) + "\n"
    )

    det_df = run_detection(
        blast_tsv=str(blast_tsv),
        assay_table=str(assay_csv),
        fna_path=str(fna_path),
        max_primer_mismatches=2,
        prime3_exact_nt=3,
        max_probe_mismatches=1,
        max_amplicon_size=500,
        store_amplicon_sequences=False,
    )

    assert len(det_df) == 1
    row = det_df.iloc[0]
    assert row['assay'] == 'VhPath'
    assert row['detection_call'] == 'Detected'
    assert row['n_amplicons'] == 1
    assert not row['multi_amplicon_flag']
    assert '201' in str(row['amplicon_sizes'])  # 300 - 100 + 1
    assert '100' in str(row['amplicon_starts']) and '300' in str(row['amplicon_ends'])


def test_run_detection_not_detected(tmp_path):
    """No BLAST hits → Not Detected for all assays."""
    from detect import run_detection

    assay_csv = tmp_path / "assay_table.csv"
    assay_csv.write_text(
        "assay,probe,fwd,rev\n"
        "VhPath,ACGGGACAAAAAGGATGGCGAGTAC,AGCCGAGCGTTACCAGC,CGAACGCAATGATTCTCTGAGC\n"
    )
    fna_path = tmp_path / "GCF_000001.fna"
    fna_path.write_text(">contig1\n" + "A" * 400 + "\n")
    blast_tsv = tmp_path / "blast.tsv"
    blast_tsv.write_text("")  # empty file

    det_df = run_detection(
        blast_tsv=str(blast_tsv),
        assay_table=str(assay_csv),
        fna_path=str(fna_path),
        max_primer_mismatches=2,
        prime3_exact_nt=3,
        max_probe_mismatches=1,
        max_amplicon_size=500,
        store_amplicon_sequences=False,
    )

    assert len(det_df) == 1
    assert det_df.iloc[0]['detection_call'] == 'Not Detected'
    # merged single frame carries the joined amplicon columns (empty when no amplicons)
    for col in ('amplicon_starts', 'amplicon_ends', 'amplicon_sequences'):
        assert col in det_df.columns


# --- call_detection mixed-probe multi-amplicon ---

def test_call_mixed_probe_multi_amplicon():
    """Two amplicons, only one probe-positive → Detected, multi_flag=True."""
    amps = [
        {'probe_found': True, 'amplicon_size': 142, 'contig_id': 'c1'},
        {'probe_found': False, 'amplicon_size': 198, 'contig_id': 'c1'},
    ]
    call, n, multi, sizes, contigs = call_detection(amps)
    assert call == 'Detected'
    assert n == 2
    assert multi is True


# --- IUPAC primitives ---

def test_iupac_match_exact_bases():
    """Standard bases match themselves."""
    for b in 'ACGT':
        assert iupac_match(b, b) is True

def test_iupac_match_degenerate_valid():
    """Degenerate codes match all bases in their set."""
    assert iupac_match('R', 'A') is True
    assert iupac_match('R', 'G') is True
    assert iupac_match('Y', 'C') is True
    assert iupac_match('Y', 'T') is True
    assert iupac_match('N', 'A') is True
    assert iupac_match('N', 'T') is True

def test_iupac_match_degenerate_invalid():
    """Degenerate codes don't match bases outside their set."""
    assert iupac_match('R', 'C') is False
    assert iupac_match('R', 'T') is False
    assert iupac_match('Y', 'A') is False
    assert iupac_match('Y', 'G') is False

def test_iupac_match_gap_never_matches():
    """Gap character never matches any base."""
    assert iupac_match('A', '-') is False
    assert iupac_match('N', '-') is False

def test_count_iupac_mismatches_perfect():
    assert count_iupac_mismatches('ACGT', 'ACGT') == 0

def test_count_iupac_mismatches_degenerate_no_mismatch():
    """R opposite A or G = 0 mismatches."""
    assert count_iupac_mismatches('ACRG', 'ACAG') == 0
    assert count_iupac_mismatches('ACRG', 'ACGG') == 0

def test_count_iupac_mismatches_degenerate_mismatch():
    """R opposite C or T = 1 mismatch."""
    assert count_iupac_mismatches('ACRG', 'ACCG') == 1

def test_count_iupac_mismatches_gap_counts():
    """Gap in sseq counts as a mismatch."""
    assert count_iupac_mismatches('ACGT', 'AC-T') == 1

def test_check_3prime_exact_iupac_pass():
    """Degenerate base at 3' position matches compatible subject base."""
    assert check_3prime_exact('ACGR', 'ACGA') is True
    assert check_3prime_exact('ACGR', 'ACGG') is True

def test_check_3prime_exact_iupac_fail():
    """Degenerate base at 3' position fails on incompatible subject base."""
    assert check_3prime_exact('ACGR', 'ACGC') is False
    assert check_3prime_exact('ACGR', 'ACGT') is False

# --- filter_primer_hits with IUPAC and gapopen ---

def test_filter_primer_hits_iupac_match_accepted():
    """Degenerate primer base matching subject within threshold is accepted."""
    # Primer: AGCCGAGCGTTACCAGR (17 nt, R at end)
    # Subject: AGCCGAGCGTTACCAGA — R matches A → 0 IUPAC mismatches
    hit = make_hit(
        qend=17, mismatch=1,  # BLAST overcounts: reports 1 mismatch for R vs A
        qseq='AGCCGAGCGTTACCAGR',
        sseq='AGCCGAGCGTTACCAGA',
        sstart=100, send=116,
    )
    contigs = {'contig1': contig_embedding(hit)}
    result = filter_primer_hits(pd.DataFrame([hit]), primer_seq='AGCCGAGCGTTACCAGR',
                                 max_mismatch=0, prime3_exact=3, contigs=contigs)
    assert len(result) == 1
    assert result.iloc[0]['mismatch'] == 0  # IUPAC-corrected count


def test_filter_primer_hits_iupac_true_mismatch_rejected():
    """Degenerate primer base NOT matching subject is still a mismatch."""
    # Primer: AGCCGAGCGTTACCAGR (R at end)
    # Subject: AGCCGAGCGTTACCAGC — R does NOT match C → 1 IUPAC mismatch
    hit = make_hit(
        qend=17, mismatch=1,
        qseq='AGCCGAGCGTTACCAGR',
        sseq='AGCCGAGCGTTACCAGC',
        sstart=100, send=116,
    )
    contigs = {'contig1': contig_embedding(hit)}
    result = filter_primer_hits(pd.DataFrame([hit]), primer_seq='AGCCGAGCGTTACCAGR',
                                 max_mismatch=0, prime3_exact=3, contigs=contigs)
    assert len(result) == 0


def test_filter_primer_hits_gapopen_rejected():
    """Hits with gapopen > 0 are rejected regardless of mismatch count."""
    hit = {**make_hit(qend=17, mismatch=0,
                      qseq='AGCCGAGCGTTACCAGC',
                      sseq='AGCCGAGCGTTACCAGC',
                      sstart=100, send=116),
           'gapopen': 1}
    contigs = {'contig1': contig_embedding(hit)}
    result = filter_primer_hits(pd.DataFrame([hit]), primer_seq='AGCCGAGCGTTACCAGC',
                                 max_mismatch=2, prime3_exact=3, contigs=contigs)
    assert len(result) == 0


# --- call_detection: probe-free assays ---

def test_call_detection_no_probe_with_amplicon():
    """has_probe=False: any valid amplicon → Detected (not Primer Only)."""
    amps = [{'contig_id': 'c1', 'amplicon_size': 200, 'probe_found': False,
             'amplicon_start': 100, 'amplicon_end': 299}]
    call, n, multi, sizes, contigs = call_detection(amps, has_probe=False)
    assert call == 'Detected'
    assert n == 1
    assert not multi


def test_call_detection_no_probe_no_amplicon():
    """has_probe=False + no amplicons → Not Detected."""
    call, n, multi, sizes, contigs = call_detection([], has_probe=False)
    assert call == 'Not Detected'
    assert n == 0


def test_call_detection_no_probe_multi_amplicon():
    """has_probe=False + 2 amplicons → Detected, multi_flag True."""
    amps = [
        {'contig_id': 'c1', 'amplicon_size': 200, 'probe_found': False,
         'amplicon_start': 100, 'amplicon_end': 299},
        {'contig_id': 'c1', 'amplicon_size': 210, 'probe_found': False,
         'amplicon_start': 500, 'amplicon_end': 709},
    ]
    call, n, multi, sizes, contigs = call_detection(amps, has_probe=False)
    assert call == 'Detected'
    assert n == 2
    assert multi


def test_check_probe_empty_seq_passthrough():
    """Empty probe_seq returns amplicons unchanged (probe_found stays False)."""
    amps = [{'contig_id': 'c1', 'amplicon_start': 100, 'amplicon_end': 300,
             'amplicon_size': 201, 'fwd_mismatches': 0, 'rev_mismatches': 0,
             'probe_found': False, 'probe_mismatches': None, 'probe_strand': None,
             'amplicon_sequence': ''}]
    probe_hits = pd.DataFrame(columns=['qseqid', 'sseqid', 'pident', 'length',
                                        'mismatch', 'gapopen', 'qstart', 'qend',
                                        'sstart', 'send', 'evalue', 'bitscore',
                                        'qseq', 'sseq'])
    result = check_probe_in_amplicons(amps, probe_hits, probe_seq='',
                                      max_probe_mismatches=1, contigs={})
    assert result[0]['probe_found'] is False
    assert result[0]['probe_mismatches'] is None


# --- revcomp / load_contigs / reconstruct_full_hits ---

def test_revcomp_iupac():
    assert revcomp('ACGT') == 'ACGT'
    assert revcomp('ACGTRY') == 'RYACGT'   # R<->Y complement, then reversed
    assert revcomp('AAA-C') == 'G-TTT'     # '-' passes through


def test_load_contigs(tmp_path):
    fna = tmp_path / 'g.fna'
    fna.write_text('>c1 description text\nacgt\nACGT\n>c2\nTTTT\n')
    contigs = load_contigs(str(fna))
    assert contigs == {'c1': 'ACGTACGT', 'c2': 'TTTT'}


def test_reconstruct_plus_strand_5prime_trim():
    """groEL_Valg scenario: BLAST trims 3 bases off the 5' end (2 mismatches at
    primer positions 2-3); reconstruction recovers the full 25-nt window."""
    primer = 'GATTCGGTGAAGAAGAGATGATCTC'          # 25 nt
    window = 'GTATCGGTGAAGAAGAGATGATCTC'          # genome: mismatches at pos 2,3
    contigs = {'c1': 'C' * 99 + window + 'C' * 50}  # window at 100..124
    hit = make_hit(sseqid='c1', qstart=4, qend=25, length=22, mismatch=0,
                   sstart=103, send=124,
                   qseq=primer[3:], sseq=window[3:])
    out = reconstruct_full_hits(pd.DataFrame([hit]), 25, contigs)
    assert len(out) == 1
    row = out.iloc[0]
    assert row['sseq'] == window
    assert (row['sstart'], row['send']) == (100, 124)
    assert (row['qstart'], row['qend']) == (1, 25)
    assert count_iupac_mismatches(primer, row['sseq']) == 2


def test_reconstruct_minus_strand_3prime_trim():
    """Tall_dnaJ scenario on the minus strand: BLAST trims the last 4 query
    bases (penultimate mismatch); reconstruction recovers them from the genome."""
    primer = 'AGCAGCTTATGACCAATACGCC'             # 22 nt
    window = 'AGCAGCTTATGACCAATACGGC'             # genome (primer-oriented), pos 21 C->G
    # minus-strand: genome plus-strand carries revcomp(window) at 200..221
    contigs = {'c1': 'A' * 199 + revcomp(window) + 'A' * 30}
    # BLAST aligned query 1..18 only; primer-oriented subject coords:
    # query base 1 at subject 221 (sstart), query base 18 at subject 204 (send)
    hit = make_hit(sseqid='c1', qstart=1, qend=18, length=18, mismatch=0,
                   sstart=221, send=204,
                   qseq=primer[:18], sseq=window[:18])
    out = reconstruct_full_hits(pd.DataFrame([hit]), 22, contigs)
    assert len(out) == 1
    row = out.iloc[0]
    assert row['sseq'] == window
    assert (row['sstart'], row['send']) == (221, 200)
    assert count_iupac_mismatches(primer, row['sseq']) == 1
    assert check_3prime_exact(primer, row['sseq'], 2) is False  # penultimate mismatch


def test_reconstruct_pads_past_contig_edge():
    """Extension beyond the contig start is padded with '-', which downstream
    filters count as mismatches."""
    primer = 'AGCCGAGCGTTACCAGC'                  # 17 nt
    contigs = {'c1': primer[2:] + 'G' * 30}       # only last 15 primer bases present, at 1..15
    hit = make_hit(sseqid='c1', qstart=3, qend=17, length=15, mismatch=0,
                   sstart=1, send=15, qseq=primer[2:], sseq=primer[2:])
    out = reconstruct_full_hits(pd.DataFrame([hit]), 17, contigs)
    assert len(out) == 1
    row = out.iloc[0]
    assert row['sseq'] == '--' + primer[2:]
    assert count_iupac_mismatches(primer, row['sseq']) == 2


def test_reconstruct_dedupes_identical_windows():
    """Two partial hits from different seeds of the same site collapse to one row."""
    primer = 'AGCCGAGCGTTACCAGC'
    contigs = {'c1': 'A' * 99 + primer + 'A' * 30}
    h1 = make_hit(sseqid='c1', qstart=1, qend=12, length=12, sstart=100, send=111,
                  qseq=primer[:12], sseq=primer[:12])
    h2 = make_hit(sseqid='c1', qstart=6, qend=17, length=12, sstart=105, send=116,
                  qseq=primer[5:], sseq=primer[5:])
    out = reconstruct_full_hits(pd.DataFrame([h1, h2]), 17, contigs)
    assert len(out) == 1


def test_reconstruct_unknown_contig_dropped():
    hit = make_hit(sseqid='no_such_contig')
    out = reconstruct_full_hits(pd.DataFrame([hit]), 17, {'c1': 'ACGT'})
    assert len(out) == 0


def test_run_detection_multiple_assays_no_state_leak(tmp_path):
    """Regression: the loaded genome dict must survive across assay iterations
    (call_detection's returned contig-id string once shadowed it, breaking
    every assay after the first)."""
    from detect import run_detection

    fwd_seq = "AGCCGAGCGTTACCAGC"
    rev_seq = "CGAACGCAATGATTCTCTGAGC"
    probe_seq = "ACGGGACAAAAAGGATGGCGAGTAC"

    assay_csv = tmp_path / "assay_table.csv"
    assay_csv.write_text(
        "assay,probe,fwd,rev\n"
        f"AssayA,{probe_seq},{fwd_seq},{rev_seq}\n"
        f"AssayB,{probe_seq},{fwd_seq},{rev_seq}\n"
    )

    seq = ("T" * 99 + fwd_seq + "T" * 33 + probe_seq + "T" * 104
           + revcomp(rev_seq) + "T" * 100)
    fna_path = tmp_path / "GCF_000002.fna"
    fna_path.write_text(f">contig1\n{seq}\n")

    rows = []
    for assay in ("AssayA", "AssayB"):
        rows.append("\t".join([f"{assay}_fwd", "contig1", "17", "0", "1", "17", "100", "116"]))
        rows.append("\t".join([f"{assay}_rev", "contig1", "22", "0", "1", "22", "300", "279"]))
        rows.append("\t".join([f"{assay}_probe", "contig1", "25", "0", "1", "25", "150", "174"]))
    blast_tsv = tmp_path / "blast.tsv"
    blast_tsv.write_text("\n".join(rows) + "\n")

    det_df = run_detection(
        blast_tsv=str(blast_tsv),
        assay_table=str(assay_csv),
        fna_path=str(fna_path),
        max_primer_mismatches=2,
        prime3_exact_nt=3,
        max_probe_mismatches=1,
        max_amplicon_size=500,
        store_amplicon_sequences=False,
    )

    assert len(det_df) == 2
    assert list(det_df['detection_call']) == ['Detected', 'Detected']


# --- BLAST input contract: column set and gzip handling ---------------------

def test_blast_cols_match_the_outfmt_in_blast_smk():
    """BLAST_COLS must match -outfmt in workflow/rules/blast.smk exactly.

    A silent drift here mislabels every column of every hit, so pin it.
    """
    smk = (Path(__file__).parent.parent / "workflow" / "rules" / "blast.smk").read_text()
    outfmt = smk.split('-outfmt "6 ', 1)[1].split('"', 1)[0].split()
    assert outfmt == BLAST_COLS


def test_load_blast_results_reads_gzipped_output(tmp_path):
    import gzip
    p = tmp_path / "blast.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("\t".join(["A_fwd", "contig1", "17", "0", "1", "17", "100", "116"]) + "\n")
    df = load_blast_results(str(p))
    assert list(df.columns) == BLAST_COLS
    assert len(df) == 1
    assert df.iloc[0]['sseqid'] == 'contig1'
    assert int(df.iloc[0]['sstart']) == 100


def test_load_blast_results_handles_empty_gzip(tmp_path):
    """An assembly with no hits gzips to a ~20-byte header, not 0 bytes, so the
    size check does not catch it and the parse must."""
    import gzip
    p = tmp_path / "blast.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write("")
    assert p.stat().st_size > 0
    df = load_blast_results(str(p))
    assert df.empty
    assert list(df.columns) == BLAST_COLS


def test_run_detection_reads_gzipped_blast(tmp_path):
    """End-to-end through the gzipped path: same call as the plain-text case."""
    import gzip
    from detect import run_detection

    fwd_seq = "AGCCGAGCGTTACCAGC"
    rev_seq = "CGAACGCAATGATTCTCTGAGC"
    probe_seq = "ACGGGACAAAAAGGATGGCGAGTAC"

    assay_csv = tmp_path / "assay_table.csv"
    assay_csv.write_text("assay,probe,fwd,rev\n"
                         f"VhPath,{probe_seq},{fwd_seq},{rev_seq}\n")
    seq = ("T" * 99 + fwd_seq + "T" * 33 + probe_seq + "T" * 104
           + revcomp(rev_seq) + "T" * 100)
    fna_path = tmp_path / "GCF_000001.fna"
    fna_path.write_text(f">contig1\n{seq}\n")

    blast_gz = tmp_path / "blast.tsv.gz"
    with gzip.open(blast_gz, "wt") as fh:
        fh.write("\t".join(["VhPath_fwd", "contig1", "17", "0", "1", "17", "100", "116"]) + "\n")
        fh.write("\t".join(["VhPath_rev", "contig1", "22", "0", "1", "22", "300", "279"]) + "\n")
        fh.write("\t".join(["VhPath_probe", "contig1", "25", "0", "1", "25", "150", "174"]) + "\n")

    det_df = run_detection(
        blast_tsv=str(blast_gz), assay_table=str(assay_csv), fna_path=str(fna_path),
        max_primer_mismatches=2, prime3_exact_nt=3, max_probe_mismatches=1,
        max_amplicon_size=500, store_amplicon_sequences=False,
    )
    assert det_df.iloc[0]['detection_call'] == 'Detected'


def test_blast_smk_sets_max_target_seqs():
    """blastn defaults -max_target_seqs to 500; with per-assembly databases the
    subjects are contigs, and fragmented drafts exceed that. Pin that the rule
    sets it explicitly (Shah et al. 2019, doi:10.1093/bioinformatics/bty833)."""
    smk = (Path(__file__).parent.parent / "workflow" / "rules" / "blast.smk").read_text()
    assert "-max_target_seqs" in smk
    assert "blast_max_target_seqs" in smk


def test_config_template_uses_permissive_blast_seeding():
    """word_size 7 silently hides mismatched short oligos (a 13-mer probe with
    one mismatch is unseedable 7.7% of the time); word_size 4 does not. Pin the
    shipped default so a future 'optimisation' cannot reintroduce the blind spot
    without failing a test."""
    import yaml
    cfg = yaml.safe_load(
        (Path(__file__).parent.parent / "config" / "config.yaml").read_text()
    )
    assert cfg["blast_word_size"] == 4
    assert cfg["blast_max_target_seqs"] >= 50000


def test_word_size_4_seeds_every_realistic_mismatch_placement():
    """The property that motivates word_size 4: a hit needs one exact run of at
    least word_size. Verify no placement of up to 2 mismatches in a 13-nt or
    longer oligo can break every run below 4 — and that word_size 7 does not
    have that property."""
    from itertools import combinations

    def worst_case_longest_run(oligo_len, n_mismatch):
        best = oligo_len
        for pos in combinations(range(oligo_len), n_mismatch):
            runs, prev = [], -1
            for p_ in list(pos) + [oligo_len]:
                runs.append(p_ - prev - 1)
                prev = p_
            best = min(best, max(runs))
        return best

    for oligo_len in (13, 17, 20, 25):
        for n_mismatch in (1, 2):
            assert worst_case_longest_run(oligo_len, n_mismatch) >= 4
    # the blind spot word_size 7 leaves, made explicit
    assert worst_case_longest_run(13, 1) < 7
    assert worst_case_longest_run(17, 2) < 7


def test_run_detection_excel_table_mac_roman_and_trailing_space(tmp_path):
    """An Excel-for-Mac table (Mac Roman, CRLF) with a trailing space after a
    primer must still detect: the space used to be counted as a primer base and
    fail the 3'-exact check, silently producing Not Detected."""
    from detect import run_detection

    fwd_seq = "AGCCGAGCGTTACCAGC"
    rev_seq = "CGAACGCAATGATTCTCTGAGC"
    assay_csv = tmp_path / "assay_table.csv"
    assay_csv.write_bytes(
        ("assay,probe,fwd,rev,notes\r\n"
         f"VhPath,,{fwd_seq} ,{rev_seq},60 °C – Fernández\r\n")
        .encode("mac_roman"))
    seq = "T" * 99 + fwd_seq + "T" * 162 + revcomp(rev_seq) + "T" * 100
    fna_path = tmp_path / "GCF_000001.fna"
    fna_path.write_text(f">contig1\n{seq}\n")
    blast_tsv = tmp_path / "blast.tsv"
    blast_tsv.write_text(
        "\t".join(["VhPath_fwd", "contig1", "17", "0", "1", "17", "100", "116"]) + "\n" +
        "\t".join(["VhPath_rev", "contig1", "22", "0", "1", "22", "300", "279"]) + "\n")
    with pytest.warns(UserWarning):
        det_df = run_detection(
            blast_tsv=str(blast_tsv), assay_table=str(assay_csv), fna_path=str(fna_path),
            max_primer_mismatches=0, prime3_exact_nt=3, max_probe_mismatches=1,
            max_amplicon_size=500, store_amplicon_sequences=False)
    assert det_df.iloc[0]['detection_call'] == 'Detected'
