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
from . import thermo as T
from . import nn as NN

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


def _sites(oligo, seq, clamp, ext=6, max_mm=5, clamp_n=2):
    """Plus-strand sites where `oligo` could PRIME: the 3'-terminal `ext` nt match exactly (the
    extension-critical zone a polymerase needs to start synthesis) AND the full-length oligo has
    <= max_mm mismatches (a coarse annealing-stability budget). Mismatches further from the 3' end
    are tolerated; a mismatch within the last `ext` nt blocks priming, so that site is NOT returned.
    (`ext` is the extension gate, not a homology seed — a primer that matches at 95% but mismatches
    2 nt from its 3' end will not appear here; use _best_hit to diagnose WHY a primer failed.)
    clamp='right': 3' end is the right terminus (forward primer). clamp='left': 3' end is the left
    terminus (reverse primer, given as rc(R)). Returns [{start, mm, ident, q3}]."""
    L = len(oligo)
    if L == 0:
        return []
    sl = min(ext, L)
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
        if len(out) >= 4000:                       # pathological low-complexity template: cap collection
            break
    out.sort(key=lambda s: (s["mm"], s["start"]))  # best-by-mismatch first, so a downstream cap keeps the best
    return out


def _best_hit(oligo, seq, clamp, ext=6):
    """Permissive best ungapped placement of the FULL oligo anywhere on the plus strand, allowing
    mismatches at ANY position including the 3' end. Used only to explain a no-product result: is a
    primer's region present-but-mismatched, or truly absent from this record? Anchors on up to three
    spread k-mers (5' / middle / 3') so one or two mismatches can't hide an otherwise-homologous site.
    Returns {ident, mm, dist3p, clean3p, win} or None (no homologous region found).
      dist3p  = nt distance of the closest mismatch to the priming 3' end (None if a perfect match)
      clean3p = True if no mismatch falls within the 3'-terminal `ext` nt (i.e. it would extend)
      win     = the aligned plus-strand template window (oligo length)."""
    L = len(oligo)
    if L < ext or len(seq) < L:
        return None
    k = min(10, L)
    offs = sorted(set([0, max(0, (L - k) // 2), L - k]))
    starts = set()
    for off in offs:
        try:
            for m in re.finditer("(?=(" + _regex(oligo[off:off + k]) + "))", seq):
                st = m.start() - off
                if 0 <= st and st + L <= len(seq):
                    starts.add(st)
                    if len(starts) >= 6000:
                        break
        except re.error:
            continue
        if len(starts) >= 6000:
            break
    if not starts:
        return None
    best = None
    for st in starts:
        win = seq[st:st + L]
        mm = L - _matches(oligo, win)
        if best is None or mm < best[0]:
            best = (mm, win)
    mm, win = best
    mism = [i for i in range(L) if not _match(oligo[i], win[i])]
    dists = [(L - 1 - i) if clamp == "right" else i for i in mism]   # distance from the priming 3' end
    dist3p = min(dists) if dists else None
    return dict(ident=round(100.0 * (L - mm) / L, 1), mm=mm, dist3p=dist3p,
                clean3p=(dist3p is None or dist3p >= ext), win=win)


def _best_ident(sites):
    return max((s["ident"] for s in sites), default=0.0)


def _probe_scan(probe, region, min_ident):
    P = probe.upper().strip()
    L = len(P)
    if not L or len(region) < L:
        return 0.0, False, None
    best = 0.0; best_win = None
    for oligo, orient in ((P, "+"), (_rc_iupac(P), "-")):
        for i in range(0, len(region) - L + 1):
            ident = 100.0 * _matches(oligo, region[i:i + L]) / L
            if ident > best:
                best = ident
                best_win = region[i:i + L] if orient == "+" else _rc_iupac(region[i:i + L])  # in P's 5'->3' frame
                if best >= 100.0:
                    break
        if best >= 100.0:
            break
    return round(best, 1), best >= min_ident, best_win


def amplify(forward, reverse, probe="", seq="", max_mm=5, clamp_n=2,
            min_product=40, max_product=3000, min_probe_ident=85.0, require_3prime=True):
    """In-silico PCR of one (F, R, probe) against one template. Returns a result dict."""
    seq = re.sub(r"[^ACGTN]", "N", (seq or "").upper().replace("U", "T"))
    _ok = r"[^ACGTRYSWKMBDHVN]"                    # restrict primers/probe to the IUPAC alphabet so a stray
    # strip_mods first: an LNA/IDT-order oligo (e.g. '/56-FAM/CTTA+CA+A+GATAT+CC+ACCACA/3IABkFQ/') must collapse
    # to its bare bases. Otherwise the '+' / mod codes would each become 'N' below, frame-shifting the probe
    # and producing a spurious uniform "probe weak" identity on every isolate.
    F = re.sub(_ok, "N", T.strip_mods(forward or "").upper().replace("U", "T").strip())   # regex metachar (| . ( etc.)
    R = re.sub(_ok, "N", T.strip_mods(reverse or "").upper().replace("U", "T").strip())   # can't corrupt the seed pattern
    _probe0 = probe or ""                                  # keep original (with '+' LNA / IUPAC degeneracy) for thermo
    probe = re.sub(_ok, "N", T.strip_mods(probe or "").upper().replace("U", "T").strip())
    if not F or not R or not seq:
        return dict(amplifies=False, product=None, f_ident=0.0, r_ident=0.0,
                    f_mm=None, r_mm=None, probe_ident=None, probe_binds=False, n_products=0,
                    f_win=None, r_win=None, p_win=None)
    rcR = _rc_iupac(R)
    # cap sites per primer so a low-complexity / degenerate primer with thousands of hits cannot make the
    # F x R product search blow up; _sites returns best-by-mismatch first, so the cap keeps the best matches
    fs = _sites(F, seq, "right", max_mm=max_mm, clamp_n=clamp_n)[:600]
    rs = _sites(rcR, seq, "left", max_mm=max_mm, clamp_n=clamp_n)[:600]
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
    # binding windows in each oligo's own 5'->3' frame, so the caller can read per-position template variation
    # (F as-is on the plus strand; R is the reverse-complement of the plus-strand region rc(R) annealed to)
    def _fwin(site):
        return seq[site["start"]: site["start"] + len(F)] if site else None
    def _rwin(site):
        return _rc_iupac(seq[site["start"]: site["start"] + len(R)]) if site else None
    def _therm(oligo, win):
        # predicted oligo:template duplex stability (Tm, ΔTm vs perfect, fraction bound at anneal),
        # mismatch/LNA/degeneracy-aware. mismatch_params parses '+'/mods/IUPAC and checks base-count vs
        # the window itself; None for ambiguous (N) windows or unscoreable oligos.
        if not oligo or not win:
            return None
        try:
            return NN.mismatch_params(oligo, win)
        except Exception:
            return None
    if best is None:
        # No extensible product. Use a permissive search to report each primer's REAL best homology
        # (so a 3'-mismatched primer doesn't read as a misleading 0%) and to classify WHY there's no
        # product: absent -> no homologous region in this record (often a partial deposit / wrong locus);
        # 3prime -> region present but a mismatch in the extension zone blocks priming; else both primers
        # would individually prime but the product fell outside the size window / wrong orientation.
        fb = _best_hit(F, seq, "right"); rb = _best_hit(rcR, seq, "left")
        def _why(b):
            if b is None or b["ident"] < 70.0:
                return "absent"
            return "ok" if b["clean3p"] else "3prime"
        parts = []
        for nm, b in (("forward", fb), ("reverse", rb)):
            w = _why(b)
            if w == "absent":
                parts.append(nm + " region absent" + ((" (best %g%%)" % b["ident"]) if b else ""))
            elif w == "3prime":
                parts.append("%s 3\u2032 mismatch %d nt from end (best %g%%)" % (nm, b["dist3p"], b["ident"]))
        reason = "; ".join(parts) if parts else "primers bind but no convergent product in the size window"
        _fw = (_fwin(fs[0]) if fs else (fb["win"] if fb else None))
        _rw = (_rwin(rs[0]) if rs else (_rc_iupac(rb["win"]) if rb else None))
        return dict(amplifies=False, product=None,
                    f_ident=(fb["ident"] if fb else 0.0), r_ident=(rb["ident"] if rb else 0.0),
                    f_mm=(fb["mm"] if fb else None), r_mm=(rb["mm"] if rb else None),
                    probe_ident=None, probe_binds=False, n_products=0, reason=reason,
                    f_win=_fw, r_win=_rw, p_win=None,
                    f_therm=_therm(forward, _fw), r_therm=_therm(reverse, _rw), p_therm=None)
    _, f, r, size, lo, hi = best
    region = seq[max(0, lo - 4):hi + 4]
    amp = seq[lo:hi]                                   # the actual amplicon (between the primer 5' ends)
    p_ident, p_binds, p_win = (_probe_scan(probe, region, min_probe_ident) if probe else (None, None, None))
    amp_tm = amp_gc = None
    if "N" not in amp and 0 < len(amp) <= 20000:       # skip if the amplicon spans an ambiguous run
        try:
            amp_tm = round(T.amplicon_tm(amp), 1); amp_gc = round(T.gc_percent(amp), 1)
        except Exception:
            amp_tm = amp_gc = None
    return dict(amplifies=True, product=size,
                f_ident=f["ident"], r_ident=r["ident"], f_mm=f["mm"], r_mm=r["mm"],
                probe_ident=p_ident, probe_binds=p_binds, n_products=n_products,
                f_win=_fwin(f), r_win=_rwin(r), p_win=p_win, amp_tm=amp_tm, amp_gc=amp_gc,
                f_therm=_therm(forward, _fwin(f)), r_therm=_therm(reverse, _rwin(r)), p_therm=_therm(_probe0, p_win))
