"""Build conservative vendor order files.

The order exporter validates sequences before producing CSV/FASTA. It intentionally rejects
opaque slash-delimited modification codes in the sequence field; probe terminal labels are added
from the selected ``kind`` so a malformed value cannot silently become a different ordered oligo.
"""
import csv
import io
import re

from . import thermo as T

SCALE = {"primer": "25nm", "probe_zen": "100nm", "probe_lna": "100nm"}
PURIF = {"primer": "STD", "probe_zen": "HPLC", "probe_lna": "HPLC"}
_ALLOWED_KINDS = frozenset(SCALE)
_BASES = frozenset("ACGTRYSWKMBDHVN")


def _safe_cell(value):
    """Prevent spreadsheet applications from interpreting a user field as a formula."""
    s = str(value if value is not None else "")
    return "'" + s if s[:1] in ("=", "+", "-", "@") else s


def _clean_name(name, fallback="oligo"):
    s = " ".join(str(name or "").replace("\r", " ").replace("\n", " ").split())
    return _safe_cell(s or fallback)


def _validate_oligo(seq, kind):
    raw = re.sub(r"\s+", "", str(seq or "").upper())
    if not raw:
        raise ValueError("empty oligo sequence")
    if "/" in raw:
        raise ValueError("slash-delimited modification codes are not accepted in the sequence field")
    if len(raw) > 180:  # includes + and * notation; the stripped backbone is checked below
        raise ValueError("oligo notation is unexpectedly long")
    i = 0
    backbone = []
    while i < len(raw):
        ch = raw[i]
        if ch == "+":
            if kind != "probe_lna":
                raise ValueError("+N LNA notation is only valid for probe_lna entries")
            if i + 1 >= len(raw) or raw[i + 1] not in "ACGT":
                raise ValueError("each '+' must immediately precede one A/C/G/T LNA base")
            backbone.append(raw[i + 1]); i += 2; continue
        if ch == "*":
            if not backbone or i + 1 >= len(raw) or raw[i + 1] in "+*":
                raise ValueError("'*' must be a linkage marker between nucleotide bases")
            i += 1; continue
        if ch not in _BASES:
            raise ValueError("invalid oligo character %r" % ch)
        backbone.append(ch); i += 1
    clean, _notes, err = T.clean_seq("".join(backbone))
    if err:
        raise ValueError(err)
    if not (6 <= len(clean) <= T.MAX_OLIGO_LEN):
        raise ValueError("oligo backbone must be 6-%d nt" % T.MAX_OLIGO_LEN)
    return raw


def _wrap(kind, seq):
    if kind in ("probe_zen", "probe_lna"):
        return "/56-FAM/%s/3IABkFQ/" % seq
    return seq


def oligo_csv(entries):
    """Return IDT-style ``Name,Sequence,Scale,Purification`` CSV after strict validation."""
    if len(entries or []) > 1000:
        raise ValueError("order export is capped at 1000 oligos")
    buf = io.StringIO(newline="")
    w = csv.writer(buf)
    w.writerow(["Name", "Sequence", "Scale", "Purification"])
    for idx, e in enumerate(entries or [], 1):
        kind = str(e.get("kind", "primer"))
        if kind not in _ALLOWED_KINDS:
            raise ValueError("entry %d has unsupported kind %r" % (idx, kind))
        seq = _validate_oligo(e.get("seq"), kind)
        w.writerow([_clean_name(e.get("name"), "oligo_%d" % idx), _wrap(kind, seq),
                    SCALE[kind], PURIF[kind]])
    return buf.getvalue()


def gblock_fasta(blocks):
    """Return FASTA for unambiguous synthetic DNA fragments; unresolved bases are rejected."""
    if len(blocks or []) > 500:
        raise ValueError("gBlock export is capped at 500 fragments")
    out = []
    for idx, b in enumerate(blocks or [], 1):
        seq = re.sub(r"\s+", "", str(b.get("seq") or "").upper())
        if not seq:
            raise ValueError("gBlock %d has an empty sequence" % idx)
        if not re.fullmatch(r"[ACGT]+", seq):
            raise ValueError("gBlock %d contains unresolved or invalid bases; resolve to A/C/G/T before ordering" % idx)
        name = re.sub(r"[^A-Za-z0-9_.+-]+", "_", str(b.get("name") or "gblock_%d" % idx)).strip("_")
        out.append(">%s\n%s" % (name or "gblock_%d" % idx, seq))
    return "\n".join(out)
