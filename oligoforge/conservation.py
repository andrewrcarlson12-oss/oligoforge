"""Per-base conservation of oligos across a target set, and discrimination vs an off-target set.

Oligo-anchored: each oligo is placed at its best ungapped offset on each sequence
(both strands), and the matched residues are tallied per oligo position. No full
MSA needed, and indels under a short oligo just show as a low-scoring placement.
Degenerate oligo bases (IUPAC) are honored. Positions are reported in the oligo's
own 5'->3' orientation, so the last position is the 3' end that governs priming.
"""
import statistics as _st
try:
    import numpy as _np
except Exception:
    _np = None
# sequence base -> code for the numpy fast path: A,C,G,T/U -> 0..3, N -> 4, anything else -> 5
_T = bytearray([5] * 256)
for _c, _v in ((65, 0), (67, 1), (71, 2), (84, 3), (85, 3), (78, 4)):
    _T[_c] = _v
_T = bytes(_T)

IUPAC = {
    "A": set("A"), "C": set("C"), "G": set("G"), "T": set("T"),
    "R": set("AG"), "Y": set("CT"), "S": set("GC"), "W": set("AT"),
    "K": set("GT"), "M": set("AC"), "B": set("CGT"), "D": set("AGT"),
    "H": set("ACT"), "V": set("ACG"), "N": set("ACGT"),
}
_COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


def _rc(s):
    return s.translate(_COMP)[::-1]


def _match(o, b):
    # N in a reference read is an unknown base (low-quality call), not a definite
    # mismatch — scoring it as a mismatch would unfairly penalize conservation and
    # over-state off-target discrimination.
    return b == "N" or b in IUPAC.get(o, set())


def _best_placement_py(o, seq):
    L = len(o)
    best = None
    for strand, t in (("+", seq), ("-", _rc(seq))):
        if len(t) < L:
            continue
        for i in range(0, len(t) - L + 1):
            w = t[i:i + L]
            s = sum(1 for a, b in zip(o, w) if _match(a, b))
            if best is None or s > best["score"]:
                best = dict(score=s, ident=s / L, offset=i, strand=strand, window=w)
    return best


def _allowed(oligo):
    # (L x 6) bool: does a sequence base of each code (A,C,G,T,N,other) match this oligo position?
    M = _np.zeros((len(oligo), 6), dtype=bool)
    for k, o in enumerate(oligo):
        s = IUPAC.get(o, set())
        for ci, b in enumerate("ACGT"):
            if b in s:
                M[k, ci] = True
        M[k, 4] = True   # N in the sequence is an unknown base -> never penalise (matches _match)
    return M


def _best_strand(codes, M, L):
    n = len(codes) - L + 1
    if n <= 0:
        return None
    score = _np.zeros(n, dtype=_np.int32)
    for k in range(L):
        score += M[k][codes[k:k + n]]
    i = int(_np.argmax(score))
    return i, int(score[i])


def best_placement(oligo, seq):
    """Best ungapped placement of oligo on seq over both strands.
    Returns dict(score, ident, offset, strand, window) or None if seq too short.
    numpy fast path gives the identical result to the reference scan but in C — a long
    mitogenome reference would otherwise make the pure-Python scan O(len x L) per oligo,
    which is what timed out (HTTP 502) on whole-mitogenome targets like cox1/cox3."""
    o = oligo.upper()
    su = seq.upper().replace("U", "T")
    L = len(o)
    if _np is None:
        return _best_placement_py(o, su)
    M = _allowed(o)
    out = None
    for strand, t in (("+", su), ("-", _rc(su))):
        if len(t) < L:
            continue
        codes = _np.frombuffer(t.encode("ascii", "replace").translate(_T), dtype=_np.uint8)
        r = _best_strand(codes, M, L)
        if r is None:
            continue
        i, sc = r
        if out is None or sc > out[1]:
            out = (strand, sc, i, t)
    if out is None:
        return None
    strand, sc, i, t = out
    return dict(score=sc, ident=sc / L, offset=i, strand=strand, window=t[i:i + L])


def conservation(oligo, sequences, min_ident=0.6):
    """Per-position match fraction of the oligo across sequences that contain the region."""
    o = oligo.upper()
    L = len(o)
    counts = [dict() for _ in range(L)]
    idents = []
    placed = 0
    for seq in sequences:
        p = best_placement(o, seq)
        if not p or p["ident"] < min_ident:
            continue
        placed += 1
        idents.append(p["ident"])
        for k, b in enumerate(p["window"]):
            counts[k][b] = counts[k].get(b, 0) + 1
    per_pos = []
    for k in range(L):
        tot = sum(counts[k].values()) or 1
        match = sum(v for b, v in counts[k].items() if _match(o[k], b))
        major = max(counts[k].items(), key=lambda x: x[1]) if counts[k] else ("-", 0)
        per_pos.append(dict(pos=k + 1, oligo=o[k], pct_match=round(100 * match / tot, 1),
                            major=major[0], dist={b: v for b, v in sorted(counts[k].items())}))
    last5 = per_pos[-5:] if L >= 5 else per_pos
    return dict(n_placed=placed, n_input=len(sequences),
                mean_ident=round(100 * _st.mean(idents), 1) if idents else 0.0,
                min_pct_match=min((p["pct_match"] for p in per_pos), default=0.0),
                worst_3prime=min((p["pct_match"] for p in last5), default=0.0),
                per_pos=per_pos)


def _gapped_ident(oligo, window):
    """Best LOCAL gapped identity (fraction over the oligo) of the oligo against a
    short subject window, plus the gap count. Used as a cross-check on the ungapped
    placement: an off-target that differs only by a small in-footprint indel matches
    closely once a gap is allowed, which the ungapped scan over-counts as mismatches
    and so over-states discrimination. Returns (pct_ident, n_gaps)."""
    o = oligo.upper().replace("U", "T")
    w = (window or "").upper().replace("U", "T")
    if not o or not w:
        return 0.0, 0
    try:
        from Bio.Align import PairwiseAligner
        al = PairwiseAligner(mode="local")
        al.match_score = 2; al.mismatch_score = -1
        al.open_gap_score = -3; al.extend_gap_score = -1
        a = al.align(w, o)[0]
        c = a.counts()
        return round(100 * c.identities / len(o), 1), c.gaps
    except Exception:
        return 0.0, 0


def discrimination(oligo, off_sequences):
    """How well the oligo MISmatches an off-target (homologous) set.
    Higher min_mismatch / lower identity = cleaner discrimination."""
    o = oligo.upper()
    L = len(o)
    rows = []
    for seq in off_sequences:
        p = best_placement(o, seq)
        if not p:
            continue
        mm = [k + 1 for k, (a, b) in enumerate(zip(o, p["window"])) if not _match(a, b)]
        mm3 = [m for m in mm if m > L - 5]
        tstr = seq.upper().replace("U", "T")
        tstr = tstr if p["strand"] == "+" else _rc(tstr)
        win = tstr[max(0, p["offset"] - 3): p["offset"] + L + 3]
        gi, gg = _gapped_ident(o, win)
        rows.append(dict(ident=round(100 * p["ident"], 1), n_mismatch=len(mm),
                         mismatch_pos=mm, n_3prime_mismatch=len(mm3),
                         gapped_ident=gi, gapped_gaps=gg))
    if not rows:
        return dict(n=0, note="oligo not located in any off-target sequence")
    offmis = [0] * L
    for r in rows:
        for m in r["mismatch_pos"]:        # m is 1-based oligo position
            offmis[m - 1] += 1
    off_mismatch_frac = [round(c / len(rows), 3) for c in offmis]
    masked = [r for r in rows if r["gapped_gaps"] > 0 and r["gapped_ident"] >= 85.0
              and r["gapped_ident"] - r["ident"] >= 10.0]
    return dict(n=len(rows),
                median_ident=round(_st.median(r["ident"] for r in rows), 1),
                max_ident=max(r["ident"] for r in rows),
                min_mismatch=min(r["n_mismatch"] for r in rows),
                min_3prime_mismatch=min(r["n_3prime_mismatch"] for r in rows),
                off_mismatch_frac=off_mismatch_frac,
                indel_masked=len(masked),
                gapped_max_ident=max((r["gapped_ident"] for r in rows), default=0.0),
                rows=rows[:25])


def analyze(oligos, target_seqs, off_seqs=None, min_ident=0.6):
    """oligos: {label: seq}. Runs conservation over targets and discrimination over off-targets."""
    out = {}
    for label, seq in oligos.items():
        rec = dict(conservation=conservation(seq, target_seqs, min_ident))
        if off_seqs:
            rec["discrimination"] = discrimination(seq, off_seqs)
        out[label] = rec
    return out
