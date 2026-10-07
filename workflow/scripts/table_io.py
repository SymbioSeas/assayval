"""Robust readers for the user-edited input tables (assay_table.csv, metadata.csv).

These files are typically edited in Excel, which saves plain "CSV" in the
platform's legacy codepage — Mac Roman on macOS, Windows-1252 on Windows — and,
in many European locales, with ';' as the delimiter. A strict UTF-8 reader then
crashes on the first en dash or accented author name in a free-text column.
Every reader of these tables goes through this module so the pipeline accepts
whatever Excel produced, warns when it had to guess, and validates the assay
table up front instead of failing silently later.

Imported by the workflow scripts, which run with this directory on sys.path.
"""
import csv
import io
import re
import sys
import warnings
from pathlib import Path

from oligo import parse_oligo, OligoError

ASSAY_REQUIRED_COLS = ('assay', 'fwd', 'rev', 'probe')
# Optional per-assay overrides of the config-wide detection thresholds. A blank
# cell (or absent column) means "use the config value".
THRESHOLD_COLS = ('max_primer_mismatches', 'prime3_exact_nt',
                  'max_probe_mismatches', 'max_amplicon_size')

# Typographic characters that legitimately appear in free text. Used to judge
# which legacy codepage decodes a file into plausible text.
_TYPOGRAPHIC = set('–—‘’‚“”„…•'
                   '·°±×÷µ§¶©®'
                   '™€£¥¢«»‹› ')

_SAVE_HINT = ('Re-save it from Excel as "CSV UTF-8 (Comma delimited)" to silence '
              'this warning.')


def _plausibility(text: str) -> int:
    """Score how much a decoding's non-ASCII characters look like real text.

    The two candidate codepages map the same high bytes to different
    characters: Mac Roman's 0xD0 is an en dash where Windows-1252 has 'Ð', and
    Windows-1252's 0xE9 'é' is Mac Roman's 'È'. Real text puts letters inside
    words (lowercase after lowercase) and punctuation between words, so each
    non-ASCII character scores +1 when it fits that pattern and -1 otherwise.
    """
    score = 0
    n = len(text)
    for i, c in enumerate(text):
        if ord(c) < 128:
            continue
        prev = text[i - 1] if i > 0 else ' '
        nxt = text[i + 1] if i + 1 < n else ' '
        in_word = prev.isalpha() and nxt.isalpha()
        touches_letter = prev.isalpha() or nxt.isalpha()
        if c.isalpha():
            if c.isupper() and prev.isalpha() and prev.islower():
                score -= 1          # 'JosÈ' — uppercase mid-word
            elif touches_letter:
                score += 1
            else:
                score -= 1          # '14Ð15' — a lone letter between digits
        elif c in _TYPOGRAPHIC:
            score += -1 if in_word else 1   # 'Caama–o' vs '14–15'
        else:
            score -= 1              # '‡', '¡', '¸' ... rarely intended
    return score


def decode_table_bytes(data: bytes, source: str = 'table') -> tuple[str, str]:
    """Decode a text table, returning (text, encoding name).

    Order: UTF-8 (BOM stripped) → UTF-16 (BOM required) → the better-scoring
    of Windows-1252 / Mac Roman → latin-1 (cannot fail). Any non-Unicode
    result emits a UserWarning naming the guess.
    """
    if data.startswith((b'\xff\xfe', b'\xfe\xff')):
        return data.decode('utf-16'), 'utf-16'
    try:
        return data.decode('utf-8-sig'), 'utf-8'
    except UnicodeDecodeError:
        pass
    candidates = []
    # Prefer the local platform's Excel codepage when the scores tie.
    order = ['mac_roman', 'cp1252'] if sys.platform == 'darwin' else ['cp1252', 'mac_roman']
    for enc in order:
        try:
            text = data.decode(enc)
        except UnicodeDecodeError:  # cp1252 leaves 5 bytes undefined
            continue
        candidates.append((_plausibility(text), -len(candidates), text, enc))
    if candidates:
        _, _, text, enc = max(candidates)
    else:
        text, enc = data.decode('latin-1'), 'latin-1'
    label = {'mac_roman': 'Mac Roman (Excel for Mac "CSV")',
             'cp1252': 'Windows-1252 (Excel for Windows "CSV")'}.get(enc, enc)
    warnings.warn(f"{source} is not UTF-8; read it as {label} [{enc}]. {_SAVE_HINT}",
                  UserWarning, stacklevel=2)
    return text, enc


def read_table_text(path) -> str:
    path = Path(path)
    text, _ = decode_table_bytes(path.read_bytes(), source=path.name)
    return text


def sniff_delimiter(text: str) -> str:
    """',' unless the header row is clearly ';'- or tab-delimited."""
    header = text.splitlines()[0] if text else ''
    counts = {d: header.count(d) for d in (',', ';', '\t')}
    best = max(counts, key=counts.get)
    return best if counts[best] > 0 else ','


def read_csv_rows(path) -> list[dict]:
    """Read a delimited table into dicts with whitespace-trimmed keys/values.

    Trimming uses str.strip(), which also removes the non-breaking spaces Excel
    sometimes leaves in cells. Fully blank rows (Excel's trailing ',,,') are
    dropped; short rows are padded with ''.
    """
    text = read_table_text(path)
    reader = csv.reader(io.StringIO(text, newline=''), delimiter=sniff_delimiter(text))
    rows = list(reader)
    if not rows:
        return []
    header = [h.strip() for h in rows[0]]
    out = []
    for raw in rows[1:]:
        cells = [c.strip() for c in raw]
        if not any(cells):
            continue
        cells += [''] * (len(header) - len(cells))
        out.append(dict(zip(header, cells)))
    return out


def load_assay_table(path) -> list[dict]:
    """Read and validate assay_table.csv; raise ValueError naming the problem.

    Validation catches mistakes that otherwise fail silently downstream: an
    assay name containing whitespace is truncated by BLAST at the space (so its
    hits never match the assay), and any non-IUPAC character in an oligo — a
    stray dash, space or smart quote — is counted as a mismatching base.
    """
    rows = read_csv_rows(path)
    name = Path(path).name
    if not rows:
        raise ValueError(f"{name}: no assays found")
    missing = [c for c in ASSAY_REQUIRED_COLS if c not in rows[0]]
    if missing:
        raise ValueError(f"{name}: Assay table missing required columns: "
                         f"{', '.join(missing)} (found: {', '.join(rows[0])})")
    seen = set()
    for i, r in enumerate(rows, start=2):  # line 1 is the header
        assay = r['assay']
        where = f"{name} line {i}"
        if not assay:
            raise ValueError(f"{where}: empty assay name")
        if re.search(r'\s', assay):
            raise ValueError(f"{where}: assay name '{assay}' contains whitespace; "
                             "use e.g. underscores instead")
        if assay in seen:
            raise ValueError(f"{where}: duplicate assay name '{assay}'")
        seen.add(assay)
        for col in ('fwd', 'rev', 'probe'):
            try:
                bases = parse_oligo(r[col]).bases
            except OligoError as e:
                raise ValueError(f"{where}: assay '{assay}' {col}: {e}") from None
            if not bases and col != 'probe':
                raise ValueError(f"{where}: assay '{assay}' {col} sequence is empty")
        for col in THRESHOLD_COLS:
            try:
                _threshold_value(r.get(col, ''))
            except ValueError:
                raise ValueError(
                    f"{where}: assay '{assay}' {col} must be blank or a "
                    f"non-negative whole number, got {r[col]!r}") from None
    return rows


def _threshold_value(cell: str):
    """None for a blank cell, else a non-negative int ('2.0' accepted, since
    Excel may write whole numbers that way); ValueError otherwise."""
    cell = (cell or '').strip()
    if not cell:
        return None
    value = float(cell)
    if value < 0 or not value.is_integer():
        raise ValueError(cell)
    return int(value)


def assay_thresholds(row: dict, defaults: dict) -> dict:
    """Effective detection thresholds for one assay: its non-blank
    THRESHOLD_COLS cells override the config-wide `defaults`."""
    out = dict(defaults)
    for col in THRESHOLD_COLS:
        value = _threshold_value(row.get(col, ''))
        if value is not None:
            out[col] = value
    return out


def read_metadata(path):
    """Read metadata.csv into a DataFrame with the same encoding/delimiter handling."""
    import pandas as pd
    text = read_table_text(path)
    df = pd.read_csv(io.StringIO(text), sep=sniff_delimiter(text))
    df.columns = [str(c).strip() for c in df.columns]
    return df
