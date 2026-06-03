"""Isolate validation: mismatch-tolerant in-silico PCR against one isolate sequence.

This is NOT the BLAST-based specificity.in_silico_pcr (which scores a pair across all of
nt). Here the caller already has a specific isolate sequence in hand (one strain genome or
gene record) and asks a direct question: would THIS pair amplify THIS template, at what
amplicon size, with what primer/probe identity and 3'-end fidelity? Run over a panel of
target isolates it measures inclusivity (sensitivity); run over near-neighbours it measures
exclusivity (cross-reactivity). The route scans one genome at a time and frees it, so peak
memory is a single record even for a 40-isolate panel.

Method: seed-and-extend. Find exact occurrences of a 3'-anchored seed of each primer (the
priming end is what matters), then verify the full primer alignment allowing up to max_mm
mismatches and requiring the terminal clamp_n bases to match (no polymerase extension off a
3' mismatch). The reverse primer is searched on the plus strand as its reverse complement;
its priming 3' end is then the LEFT terminus of that match. A product is a convergent F/R
pair within the size window. The probe is scanned across the amplicon on both strands with
no 3' requirement (it does not prime).
"""
import re
from .specificity import _rc_iupac

# IUPAC -> the set of plain bases it matches
_DEG = {"A": "A", "C": "C", "G": "G", "T": "T",
        "R": "AG", "Y": "CT", "S": "GC", "W": "AT", "K": "GT", "M": "AC",
        "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}


def _match(o, b):
    return b in _DEG.get(o, o)


def _matches(oligo, window):
    return sum(1 for o, b in zip(oligo, window) if _match(o, b))


def _regex(seq):
    return "".join(("[" + _DEG[b] + "]" if len(_DEG.get(b, b)) > 1 else b) for b in seq)


def _sites(oligo, seq, clamp, seed_len=13, max_mm=5, clamp_n=2):
    """Plus-strand binding sites of `oligo` (5'->3' as it sits on the plus strand).
    clamp='right': priming 3' end is the right terminus (a forward primer).
    clamp='left' : priming 3' end is the left terminus (a reverse primer, given as rc(R)).
    Returns [{start, mm, ident, q3}]."""
    L = len(oligo)
    if L == 0:
        return []
    sl = min(seed_len, L)
    seed = oligo[-sl:] if clamp == "right" else oligo[:sl]
    pat = "(?=(" + _regex(seed) + "))"
    out = []
    try:
        it = re.finditer(pat, seq)
    except re.error:
        return []
    for m in it:
        s = m.start()
        start = s - (L - sl) if clamp == "right" else s
        if start < 0 or start + L > len(seq):
            continue
        win = seq[start:start + L]
        mat = _matches(oligo, win)
        mm = L - mat
        if mm > max_mm:
            continue
        if clamp == "right":
            q3 = _matches(oligo[-clamp_n:], win[-clamp_n:]) == clamp_n
        else:
            q3 = _matches(oligo[:clamp_n], win[:clamp_n]) == clamp_n
        out.append(dict(start=start, mm=mm, ident=round(100.0 * mat / L, 1), q3=q3))
    return out


def _best_ident(sites):
    return max((s["ident"] for s in sites), default=0.0)


def _probe_scan(probe, region, min_ident):
    P = probe.upper().strip()
    L = len(P)
    if not L or len(region) < L:
        return 0.0, False
    best = 0.0
    for oligo in (P, _rc_iupac(P)):
        for i in range(0, len(region) - L + 1):
            ident = 100.0 * _matches(oligo, region[i:i + L]) / L
            if ident > best:
                best = ident
                if best >= 100.0:
                    break
    return round(best, 1), best >= min_ident


def amplify(forward, reverse, probe="", seq="", max_mm=5, clamp_n=2,
            min_product=40, max_product=3000, min_probe_ident=85.0, require_3prime=True):
    """In-silico PCR of one (F, R, probe) against one template. Returns a result dict."""
    seq = re.sub(r"[^ACGTN]", "N", (seq or "").upper().replace("U", "T"))
    F = (forward or "").upper().strip()
    R = (reverse or "").upper().strip()
    if not F or not R or not seq:
        return dict(amplifies=False, product=None, f_ident=0.0, r_ident=0.0,
                    f_mm=None, r_mm=None, probe_ident=None, probe_binds=False, n_products=0)
    rcR = _rc_iupac(R)
    fs = _sites(F, seq, "right", max_mm=max_mm, clamp_n=clamp_n)
    rs = _sites(rcR, seq, "left", max_mm=max_mm, clamp_n=clamp_n)
    best = None
    n_products = 0
    for f in fs:
        if require_3prime and not f["q3"]:
            continue
        fi = f["start"]
        for r in rs:
            if require_3prime and not r["q3"]:
                continue
            rj = r["start"]
            size = (rj + len(R)) - fi
            if fi < rj and min_product <= size <= max_product:
                n_products += 1
                key = (size, f["mm"] + r["mm"])
                if best is None or key < best[0]:
                    best = (key, f, r, size, fi, rj + len(R))
    if best is None:
        # no productive amplicon — still report best primer identities so the caller can
        # see WHY a near-neighbour was rejected (e.g. "best F 78%, 3' mismatch")
        return dict(amplifies=False, product=None,
                    f_ident=_best_ident(fs), r_ident=_best_ident(rs),
                    f_mm=None, r_mm=None, probe_ident=None, probe_binds=False, n_products=0)
    _, f, r, size, lo, hi = best
    region = seq[max(0, lo - 4):hi + 4]
    p_ident, p_binds = (_probe_scan(probe, region, min_probe_ident) if probe else (None, None))
    return dict(amplifies=True, product=size,
                f_ident=f["ident"], r_ident=r["ident"], f_mm=f["mm"], r_mm=r["mm"],
                probe_ident=p_ident, probe_binds=p_binds, n_products=n_products)
