"""Per-base conservation of oligos across a target set, and discrimination vs an off-target set.

Oligo-anchored: each oligo is placed at its best ungapped offset on each sequence
(both strands), and the matched residues are tallied per oligo position. No full
MSA needed, and indels under a short oligo just show as a low-scoring placement.
Degenerate oligo bases (IUPAC) are honored. Positions are reported in the oligo's
own 5'->3' orientation, so the last position is the 3' end that governs priming.
"""
import statistics as _st

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
    return b in IUPAC.get(o, set())


def best_placement(oligo, seq):
    """Best ungapped placement of oligo on seq over both strands.
    Returns dict(score, ident, offset, strand, window) or None if seq too short."""
    o = oligo.upper()
    seq = seq.upper().replace("U", "T")
    L = len(o)
    best = None
    for strand, t in (("+", seq), ("-", _rc(seq))):
        if len(t) < L:
            continue
        for i in range(0, len(t) - L + 1):
            w = t[i:i + L]
            s = 0
            for a, b in zip(o, w):
                if _match(a, b):
                    s += 1
            if best is None or s > best["score"]:
                best = dict(score=s, ident=s / L, offset=i, strand=strand, window=w)
    return best


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
        rows.append(dict(ident=round(100 * p["ident"], 1), n_mismatch=len(mm),
                         mismatch_pos=mm, n_3prime_mismatch=len(mm3)))
    if not rows:
        return dict(n=0, note="oligo not located in any off-target sequence")
    return dict(n=len(rows),
                median_ident=round(_st.median(r["ident"] for r in rows), 1),
                max_ident=max(r["ident"] for r in rows),
                min_mismatch=min(r["n_mismatch"] for r in rows),
                min_3prime_mismatch=min(r["n_3prime_mismatch"] for r in rows),
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
