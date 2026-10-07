"""Parse a primer/probe cell from assay_table.csv into the bases AssayVal scores.

Oligos are written with vendor notation: IDT-style modification codes in
slashes (/56-FAM/, /ZEN/, /3IABkFQ/), bracketed modifications ([BHQ1]),
locked nucleic acids (+G or [G]), phosphorothioate bonds (*) and so on. Some
of these carry a base (an LNA, 5-methyl-dC, inosine, a labelled dT) and some do
not (dyes, quenchers, spacers). Deleting every token wholesale — the old
behaviour — silently removed bases and shifted the oligo, so each token is now
interpreted, and anything unrecognised is an error rather than a guess.

Single source of truth for oligo parsing: imported by table_io (validation),
prepare_oligos (BLAST query FASTA) and detect (scoring).
"""
import re
from dataclasses import dataclass

IUPAC = set('ACGTRYSWKMBDHVN')
LNA_BASES = set('ACGTU')

# IDT internal codes that stand for a base (keys lowercased).
INTERNAL_BASE = {
    'ime-dc': ('C', '5-methyl-dC'),
    'ideoxyi': ('N', 'deoxyinosine, scored as N'),
    'ideoxyu': ('T', 'deoxyuridine'),
    'ibiodt': ('T', 'biotin-dT'),
    'ifluort': ('T', 'fluorescein-dT'),
    'iammc6t': ('T', 'amino-modifier C6 dT'),
}
# IDT internal codes with no base (internal quenchers, spacers).
INTERNAL_NONBASE = {'zen', 'izen', 'tao', 'itao', 'ispc3', 'isp9', 'isp18'}

_DEFERRED = {'r': 'RNA base (rN)', 'm': "2'-O-methyl RNA base (mN)"}
_DEFERRED_HINT = ("is not supported yet: RNA and 2'-O-methyl bases (including "
                  "rhAmp/rhPCR primers) need their own detection model")


class OligoError(ValueError):
    """An oligo cell that cannot be interpreted unambiguously."""


@dataclass(frozen=True)
class Oligo:
    raw: str
    bases: str          # uppercase IUPAC: what is searched and scored
    lna: tuple          # 0-based positions in `bases` that are LNA
    notes: tuple        # what was removed/converted, for prep.log


def _deferred_error(kind: str, token: str) -> OligoError:
    return OligoError(f"{token!r} looks like a {_DEFERRED[kind]}, which "
                      f"{_DEFERRED_HINT}")


def _check_deferred_mixed_case(raw: str) -> None:
    """IDT writes RNA / 2'-O-methyl bases as rA / mA. Lowercase r and m are
    also the IUPAC codes R and M, so only flag them in a mixed-case sequence,
    where a lowercase letter among uppercase bases is clearly a prefix."""
    letters = re.sub(r'/[^/]*/|\[[^\]]*\]', '', raw)
    if not (any(c.isupper() for c in letters) and any(c.islower() for c in letters)):
        return
    m = re.search(r'([rm])([ACGU])', letters)
    if m:
        raise _deferred_error(m.group(1), m.group(0))


def parse_oligo(raw: str) -> Oligo:
    raw = (raw or '').strip()
    if raw[:1] in ('#', '=', '@'):
        raise OligoError(
            f"{raw!r} looks like an Excel formula or error value, not a sequence. "
            "Excel treats a cell starting with '+', '=', '-' or '@' as a formula; "
            "write a 5'-terminal LNA as [X] instead of +X")
    _check_deferred_mixed_case(raw)

    bases, lna, notes = [], [], []
    three_prime = None  # 3' modification seen; no base may follow it
    saw_ps = False

    def add(base, is_lna=False):
        if three_prime:
            raise OligoError(f"3' modification {three_prime} must come last "
                             "(a base follows it)")
        if is_lna:
            lna.append(len(bases))
        bases.append(base)

    i, n = 0, len(raw)
    while i < n:
        c = raw[i]
        if c == '/':
            j = raw.find('/', i + 1)
            if j == -1:
                raise OligoError(f"unclosed '/' in {raw!r}")
            code = raw[i + 1:j]
            token = f"/{code}/"
            key = code.lower()
            if key in INTERNAL_BASE:
                base, what = INTERNAL_BASE[key]
                add(base)
                notes.append(f"converted {token} -> {base} ({what})")
            elif key in INTERNAL_NONBASE:
                notes.append(f"removed internal modification {token}")
            elif code.startswith('5'):
                if bases:
                    raise OligoError(f"5' modification {token} must come first")
                notes.append(f"removed 5' modification {token}")
            elif code.startswith('3'):
                three_prime = token
                notes.append(f"removed 3' modification {token}")
            elif key.startswith('i'):
                raise OligoError(
                    f"unrecognized internal modification {token}. If it carries "
                    "no base, write it in brackets, e.g. [name]; if it stands "
                    "for a base, write that base instead")
            else:
                raise OligoError(
                    f"unrecognized modification {token}. IDT codes start with "
                    "5 (5' end), 3 (3' end) or i (internal); write other "
                    "non-base modifications in brackets, e.g. [name]")
            i = j + 1
        elif c == '[':
            j = raw.find(']', i + 1)
            if j == -1:
                raise OligoError(f"unclosed '[' in {raw!r}")
            content = raw[i + 1:j].strip()
            token = f"[{content}]"
            core = content[1:] if content.startswith('+') else content
            if len(core) == 1:
                if core.upper() not in LNA_BASES:
                    raise OligoError(
                        f"{token}: a single base in brackets marks an LNA and must "
                        "be A, C, G or T")
                add('T' if core.upper() == 'U' else core.upper(), is_lna=True)
            elif len(content) == 2 and content[0] in 'rRmM' \
                    and content[1].upper() in LNA_BASES:
                raise _deferred_error(content[0].lower(), token)
            else:
                notes.append(f"removed non-base modification {token}")
            i = j + 1
        elif c == '+':
            nxt = raw[i + 1:i + 2]
            if not nxt or nxt.upper() not in LNA_BASES:
                raise OligoError(f"'+' must be followed by A, C, G or T (LNA) "
                                 f"in {raw!r}")
            add('T' if nxt.upper() == 'U' else nxt.upper(), is_lna=True)
            i += 2
        elif c == '*':
            saw_ps = True
            i += 1
        else:
            u = c.upper()
            if u == 'I':
                add('N')
                notes.append("converted inosine I -> N")
            elif u == 'U':
                add('T')
                notes.append("converted U -> T")
            elif u in IUPAC:
                add(u)
            else:
                raise OligoError(f"invalid character {c!r} in {raw!r}; only IUPAC "
                                 "bases and modification notation are allowed")
            i += 1

    if saw_ps:
        notes.append("removed phosphorothioate bond(s) '*'")
    return Oligo(raw=raw, bases=''.join(bases), lna=tuple(lna), notes=tuple(notes))
