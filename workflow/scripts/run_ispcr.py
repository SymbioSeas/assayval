import re
import csv
import argparse
import pandas as pd
from pathlib import Path
from Bio import SeqIO


BLAST_COLS = [
    'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
    'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qseq', 'sseq',
]

DETECTION_COLS = [
    'accession', 'assay', 'detection_call', 'n_amplicons', 'multi_amplicon_flag',
    'amplicon_sizes', 'contig_ids', 'amplicon_starts', 'amplicon_ends',
    'fwd_mismatches', 'rev_mismatches', 'probe_mismatches', 'probe_strand',
    'amplicon_sequences',
]

_MOD_RE = re.compile(r'/[^/]+/|\[[^\]]+\]')


def _strip(seq: str) -> str:
    return _MOD_RE.sub('', seq)


IUPAC_BASES = {
    'A': {'A'},      'T': {'T'},      'G': {'G'},      'C': {'C'},
    'N': {'A', 'T', 'G', 'C'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
    'K': {'G', 'T'}, 'M': {'A', 'C'},
    'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'},
}


def iupac_match(q: str, s: str) -> bool:
    """True if subject base s is within the set represented by query IUPAC base q."""
    return s.upper() in IUPAC_BASES.get(q.upper(), {q.upper()})


def count_iupac_mismatches(primer_seq: str, sseq: str) -> int:
    """Count positions where primer base doesn't IUPAC-match the aligned subject base.
    Gap characters in sseq always count as mismatches.
    """
    return sum(
        1 for q, s in zip(primer_seq, sseq)
        if s == '-' or not iupac_match(q, s)
    )


def check_3prime_exact(primer_seq: str, sseq: str, n: int = 3) -> bool:
    """Return True if last n positions match by IUPAC rules (no gaps allowed)."""
    tail_q, tail_s = primer_seq[-n:], sseq[-n:]
    if '-' in tail_q or '-' in tail_s:
        return False
    return all(iupac_match(q, s) for q, s in zip(tail_q, tail_s))


COMPLEMENT = str.maketrans(
    'ACGTRYSWKMBDHVNacgtryswkmbdhvn',
    'TGCAYRSWMKVHDBNtgcayrswmkvhdbn',
)


def revcomp(seq: str) -> str:
    """Reverse-complement, IUPAC-aware; characters without a complement
    (e.g. '-') pass through unchanged."""
    return seq.translate(COMPLEMENT)[::-1]


def count_degenerate(primer_seq: str) -> int:
    """Number of IUPAC-degenerate positions in an oligo."""
    return sum(1 for b in primer_seq.upper() if b not in 'ACGT')


def load_contigs(fna_path: str) -> dict[str, str]:
    """Load all contigs of an assembly into memory, uppercased, keyed by id."""
    return {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(fna_path, 'fasta')}


def reconstruct_full_hits(hits: pd.DataFrame, primer_len: int,
                          contigs: dict[str, str]) -> pd.DataFrame:
    """Rebuild each gapless BLAST hit as a full-length oligo window.

    blastn reports local alignments, and extending an alignment through a
    terminal mismatch always lowers its score (blastn-short scores matches +1,
    mismatches -3), so an oligo whose outer bases mismatch the template comes
    back with qstart > 1 and/or qend < oligo length. Requiring full query
    coverage would silently discard such hits no matter how permissive the
    mismatch tolerance is. Instead, the un-aligned oligo ends are mapped onto
    subject coordinates and the complete window is re-extracted from the
    genome, so the IUPAC-aware mismatch count and 3'-exact filter judge every
    oligo position against real template sequence.

    Returned rows have sstart/send/sseq rewritten to the full-length window
    (sseq oriented to the query; minus-strand windows reverse-complemented),
    qstart=1 and qend=primer_len. Window positions beyond a contig edge are
    padded with '-', which downstream filters count as mismatches. Hits whose
    reconstructed windows are identical are deduplicated.

    Assumes gapless hits (filter gapopen==0 upstream): with no indels, query
    and subject coordinates map 1:1.
    """
    rows = []
    for _, h in hits.iterrows():
        contig = contigs.get(h['sseqid'])
        if contig is None:
            continue
        left_ext = int(h['qstart']) - 1
        right_ext = primer_len - int(h['qend'])
        sstart, send = int(h['sstart']), int(h['send'])
        if sstart <= send:  # plus strand
            lo = sstart - left_ext
            hi = send + right_ext
            new_sstart, new_send = lo, hi
        else:  # minus strand: query 5'->3' runs high->low on the subject
            lo = send - right_ext
            hi = sstart + left_ext
            new_sstart, new_send = hi, lo
        left_pad = max(0, 1 - lo)
        right_pad = max(0, hi - len(contig))
        window = contig[max(lo, 1) - 1:min(hi, len(contig))]
        window = '-' * left_pad + window + '-' * right_pad
        h = h.copy()
        h['sseq'] = revcomp(window) if sstart > send else window
        h['sstart'], h['send'] = new_sstart, new_send
        h['qstart'], h['qend'] = 1, primer_len
        rows.append(h)
    if not rows:
        return pd.DataFrame(columns=hits.columns)
    return pd.DataFrame(rows).drop_duplicates(subset=['sseqid', 'sstart', 'send'])


def filter_primer_hits(hits: pd.DataFrame, primer_seq: str,
                        max_mismatch: int, prime3_exact: int,
                        contigs: dict[str, str]) -> pd.DataFrame:
    """Return hits passing mismatch and 3'-exact filters over the FULL oligo.

    BLAST hits are gapless-filtered, reconstructed to full oligo length from
    the genome (see reconstruct_full_hits — local alignments cannot include
    terminal mismatches, so partial hits must be recovered rather than
    dropped), then judged by IUPAC-aware mismatch counting and the 3'-exact
    check. Strand filtering is deferred to find_valid_amplicons, which handles
    both target gene orientations (+ strand and - strand assemblies).
    """
    if hits.empty:
        return hits
    primer_len = len(primer_seq)
    hits = hits[hits['gapopen'] == 0].copy()
    if hits.empty:
        return hits
    # Cheap pre-filter: under blastn-short scoring (+1/-3) any trimmed
    # terminus of k bases contains >= k/4 BLAST-scored mismatches, so a hit
    # shorter than primer_len - 4*budget cannot pass the mismatch filter.
    # Degenerate oligo positions are BLAST-scored as mismatches yet may be
    # IUPAC matches, hence the count_degenerate allowance.
    min_len = primer_len - 4 * (max_mismatch + count_degenerate(primer_seq))
    hits = hits[hits['length'] >= min_len]
    if hits.empty:
        return hits
    hits = reconstruct_full_hits(hits, primer_len, contigs)
    if hits.empty:
        return hits
    hits['mismatch'] = hits['sseq'].apply(
        lambda s: count_iupac_mismatches(primer_seq, s)
    )
    hits = hits[hits['mismatch'] <= max_mismatch]
    if hits.empty:
        return hits
    hits = hits[hits['sseq'].apply(
        lambda s: check_3prime_exact(primer_seq, s, prime3_exact)
    )]
    return hits


def find_valid_amplicons(fwd_hits: pd.DataFrame, rev_hits: pd.DataFrame,
                          max_amplicon_size: int) -> list[dict]:
    """Find fwd+rev pairs that form a valid amplicon on the same contig.

    Valid pairs must be on opposite strands (one + strand, one - strand) and
    converging: both primers' 3' ends must fall within the amplicon span. This
    handles target genes on either the + strand (fwd +, rev -) or - strand
    (fwd -, rev +) of the assembly.
    """
    amplicons = []
    for _, fwd in fwd_hits.iterrows():
        fwd_plus = int(fwd['sstart']) < int(fwd['send'])
        for _, rev in rev_hits.iterrows():
            if fwd['sseqid'] != rev['sseqid']:
                continue
            rev_plus = int(rev['sstart']) < int(rev['send'])
            if fwd_plus == rev_plus:  # same strand → primers diverge, no amplicon
                continue
            amp_start = min(int(fwd['sstart']), int(rev['sstart']))
            amp_end = max(int(fwd['sstart']), int(rev['sstart']))
            # Convergence: both 3' ends must point inward (within amplicon bounds)
            fwd_3p = int(fwd['send'])
            rev_3p = int(rev['send'])
            if not (amp_start <= fwd_3p <= amp_end and amp_start <= rev_3p <= amp_end):
                continue
            amp_size = amp_end - amp_start + 1
            if amp_size > max_amplicon_size:
                continue
            amplicons.append({
                'contig_id': fwd['sseqid'],
                'amplicon_start': amp_start,
                'amplicon_end': amp_end,
                'amplicon_size': amp_size,
                'fwd_mismatches': int(fwd['mismatch']),
                'rev_mismatches': int(rev['mismatch']),
                'probe_found': False,
                'probe_mismatches': None,
                'probe_strand': None,
                'amplicon_sequence': '',
            })
    return amplicons


def check_probe_in_amplicons(amplicons: list[dict], probe_hits: pd.DataFrame,
                               probe_seq: str, max_probe_mismatches: int,
                               contigs: dict[str, str]) -> list[dict]:
    """For each amplicon, find a probe hit contained within it (either strand).

    Probe hits get the same treatment as primers: gapless filter, full-length
    reconstruction from the genome (recovers BLAST end-trimmed hits), then
    IUPAC-aware mismatch filtering before checking spatial containment.
    """
    if not probe_seq:
        return amplicons
    probe_len = len(probe_seq)
    valid_probe = probe_hits[probe_hits['gapopen'] == 0].copy()
    if not valid_probe.empty:
        min_len = probe_len - 4 * (max_probe_mismatches + count_degenerate(probe_seq))
        valid_probe = valid_probe[valid_probe['length'] >= min_len]
        valid_probe = reconstruct_full_hits(valid_probe, probe_len, contigs)
    if not valid_probe.empty:
        valid_probe['mismatch'] = valid_probe['sseq'].apply(
            lambda s: count_iupac_mismatches(probe_seq, s)
        )
        valid_probe = valid_probe[valid_probe['mismatch'] <= max_probe_mismatches]
    for amp in amplicons:
        contig_probe = valid_probe[valid_probe['sseqid'] == amp['contig_id']]
        for _, hit in contig_probe.iterrows():
            h_start = min(int(hit['sstart']), int(hit['send']))
            h_end = max(int(hit['sstart']), int(hit['send']))
            if h_start >= amp['amplicon_start'] and h_end <= amp['amplicon_end']:
                amp['probe_found'] = True
                amp['probe_mismatches'] = int(hit['mismatch'])
                amp['probe_strand'] = '+' if hit['sstart'] < hit['send'] else '-'
                break
    return amplicons


def call_detection(amplicons: list[dict], has_probe: bool = True) -> tuple:
    """Return (detection_call, n_amplicons, multi_flag, sizes_str, contigs_str).

    For probe-free assays (has_probe=False) any valid amplicon is a detection.
    """
    if not amplicons:
        return 'Not Detected', 0, False, '', ''
    n = len(amplicons)
    if not has_probe:
        call = 'Detected'
    else:
        call = 'Detected' if any(a['probe_found'] for a in amplicons) else 'Primer Only'
    multi = n > 1
    sizes = ';'.join(str(a['amplicon_size']) for a in amplicons)
    contigs = ';'.join(a['contig_id'] for a in amplicons)
    return call, n, multi, sizes, contigs


def extract_amplicon_sequence(contigs: dict[str, str], contig_id: str,
                               start: int, end: int) -> str:
    """Extract 1-based inclusive subsequence from the loaded contigs.

    The max(start, 1) guard protects against edge-reconstructed hit
    coordinates that fall before the contig start.
    """
    seq = contigs.get(contig_id, '')
    return seq[max(start, 1) - 1:end]


def load_blast_results(blast_tsv: str) -> pd.DataFrame:
    """Load BLAST tabular output; return empty DataFrame if file is absent or has no hits."""
    p = Path(blast_tsv)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame(columns=BLAST_COLS)
    df = pd.read_csv(p, sep='\t', header=None, names=BLAST_COLS)
    return df


def run_ispcr(blast_tsv: str, assay_table: str, fna_path: str,
              max_primer_mismatches: int, prime3_exact_nt: int,
              max_probe_mismatches: int, max_amplicon_size: int,
              store_amplicon_sequences: bool) -> pd.DataFrame:
    """
    Core detection engine. Returns one detection DataFrame with a single row per
    assay; per-amplicon positions and sequences are ';'-joined into that row.
    """
    blast = load_blast_results(blast_tsv)
    accession = Path(fna_path).stem
    # NB: named contig_seqs (not contigs) — the per-assay loop below rebinds
    # `contigs` to call_detection's joined contig-id string.
    contig_seqs = load_contigs(fna_path)

    with open(assay_table, encoding='utf-8-sig') as f:
        assays = list(csv.DictReader(f))

    detection_rows = []

    for row in assays:
        name = row['assay']
        fwd_seq = _strip(row['fwd'])
        rev_seq = _strip(row['rev'])
        probe_seq = _strip(row['probe'])
        has_probe = bool(probe_seq)

        fwd_hits_all = blast[blast['qseqid'] == f'{name}_fwd']
        rev_hits_all = blast[blast['qseqid'] == f'{name}_rev']
        probe_hits_all = blast[blast['qseqid'] == f'{name}_probe']

        fwd_hits = filter_primer_hits(fwd_hits_all, fwd_seq,
                                       max_primer_mismatches, prime3_exact_nt,
                                       contig_seqs)
        rev_hits = filter_primer_hits(rev_hits_all, rev_seq,
                                       max_primer_mismatches, prime3_exact_nt,
                                       contig_seqs)

        amplicons = find_valid_amplicons(fwd_hits, rev_hits, max_amplicon_size)
        amplicons = check_probe_in_amplicons(amplicons, probe_hits_all, probe_seq,
                                             max_probe_mismatches, contig_seqs)

        if store_amplicon_sequences:
            for amp in amplicons:
                amp['amplicon_sequence'] = extract_amplicon_sequence(
                    contig_seqs, amp['contig_id'], amp['amplicon_start'], amp['amplicon_end']
                )

        call, n, multi, sizes, contigs = call_detection(amplicons, has_probe)

        best_fwd = min((a['fwd_mismatches'] for a in amplicons), default=None)
        best_rev = min((a['rev_mismatches'] for a in amplicons), default=None)
        best_probe = min(
            (a['probe_mismatches'] for a in amplicons if a['probe_mismatches'] is not None),
            default=None,
        )
        probe_strand = next((a['probe_strand'] for a in amplicons if a['probe_strand']), None)

        starts = ';'.join(str(a['amplicon_start']) for a in amplicons)
        ends = ';'.join(str(a['amplicon_end']) for a in amplicons)
        seqs = ';'.join(a['amplicon_sequence'] for a in amplicons)

        detection_rows.append({
            'accession': accession, 'assay': name, 'detection_call': call,
            'n_amplicons': n, 'multi_amplicon_flag': multi,
            'amplicon_sizes': sizes, 'contig_ids': contigs,
            'amplicon_starts': starts, 'amplicon_ends': ends,
            'fwd_mismatches': best_fwd, 'rev_mismatches': best_rev,
            'probe_mismatches': best_probe, 'probe_strand': probe_strand,
            'amplicon_sequences': seqs,
        })

    return pd.DataFrame(detection_rows, columns=DETECTION_COLS)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--blast', required=True)
    p.add_argument('--fna', required=True)
    p.add_argument('--assay-table', required=True)
    p.add_argument('--max-primer-mismatches', type=int, default=2)
    p.add_argument('--prime3-exact-nt', type=int, default=1)
    p.add_argument('--max-probe-mismatches', type=int, default=1)
    p.add_argument('--max-amplicon-size', type=int, default=500)
    p.add_argument('--store-amplicon-sequences',
                   type=lambda x: x.lower() == 'true', default=True)
    p.add_argument('--detection-out', required=True)
    args = p.parse_args()

    det_df = run_ispcr(
        blast_tsv=args.blast, assay_table=args.assay_table, fna_path=args.fna,
        max_primer_mismatches=args.max_primer_mismatches,
        prime3_exact_nt=args.prime3_exact_nt,
        max_probe_mismatches=args.max_probe_mismatches,
        max_amplicon_size=args.max_amplicon_size,
        store_amplicon_sequences=args.store_amplicon_sequences,
    )
    det_df.to_csv(args.detection_out, index=False)


if __name__ == '__main__':
    main()
